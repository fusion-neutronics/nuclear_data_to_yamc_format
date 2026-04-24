#!/usr/bin/env python3

"""
Convert ENDF/B-VIII.0 or VII.1 to simulation-ready Arrow format.

Pipeline: download ENDF -> NJOY (Doppler broaden) -> OpenMC parse -> Arrow.

If ENDF source files already exist locally, the download step is skipped.
Search order: ./endfb-{release}-endf/ then ~/nuclear_data/endfb-{release}-endf/.
Output defaults to ~/nuclear_data/endf-b{X}.0-arrow/.
"""

import argparse
import re
import sys
from multiprocessing import Pool
from pathlib import Path

from nuclear_data_to_yamc_format import convert_neutron, convert_photon
from nuclear_data_to_yamc_format.cli import nuclide_filter, write_index
from nuclear_data_to_yamc_format.download import (
    ENDF_RELEASES, download_and_extract, find_photon_files,
)

assert sys.version_info >= (3, 9), "Python 3.9+ is required"


# Filename patterns the convert pipeline expects for each particle. Used to
# decide whether an existing local ENDF source already has the files needed
# for that particle, or whether we need to download them.
PARTICLE_GLOBS = {
    "neutron": ("n-*.endf",),
    "photon": ("photoat-*.endf", "atom-*.endf"),
}


def _has_particle_source(endf_dir, particle):
    globs = PARTICLE_GLOBS.get(particle, ())
    return bool(globs) and all(
        next(endf_dir.rglob(g), None) is not None for g in globs
    )


_NEUTRON_FILENAME_RE = re.compile(
    r"^n-\d+_(?P<sym>[A-Za-z]+)_(?P<mass>\d+)(?:m(?P<meta>\d+))?$"
)


def _nuclide_name_from_neutron_endf(path):
    """Derive the expected Arrow output name (e.g. 'Ag110_m1') from a neutron
    ENDF filename like 'n-047_Ag_110m1.endf'. Returns None if the filename
    doesn't match the expected pattern — the caller should then fall back to
    running the conversion (which parses the file for the authoritative name).
    """
    m = _NEUTRON_FILENAME_RE.match(path.stem)
    if not m:
        return None
    name = f"{m['sym']}{int(m['mass'])}"
    if m["meta"]:
        name += f"_m{m['meta']}"
    return name


def find_or_download_endf(release, particles):
    """Return the ENDF source dir, downloading any requested particle whose
    source files are missing.

    Each particle is checked independently, so a pre-existing neutron source
    does not suppress a needed photon download.
    """
    info = ENDF_RELEASES[release]
    dirname = info["dir"]

    candidates = [
        Path.cwd() / dirname,
        Path.home() / "nuclear_data" / dirname,
    ]
    endf_dir = next(
        (p for p in candidates if p.is_dir() and any(p.rglob("*.endf"))),
        Path.home() / "nuclear_data" / dirname,
    )
    download_dir = endf_dir / "_downloads"

    for particle in particles:
        if particle not in info:
            continue
        if _has_particle_source(endf_dir, particle):
            print(f"Using existing {particle} ENDF source in {endf_dir}")
            continue
        details = info[particle]
        urls = [details["base_url"] + f for f in details["files"]]
        print(f"\nDownloading {particle} data to {endf_dir}")
        download_and_extract(urls, endf_dir, download_dir)

    return endf_dir


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=type(
            "F",
            (argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter),
            {},
        ),
    )
    parser.add_argument("-d", "--destination", type=Path,
                        help="Directory to create new library in")
    parser.add_argument("-r", "--release", choices=list(ENDF_RELEASES.keys()),
                        default="viii.1", help="ENDF/B release version")
    parser.add_argument("-p", "--particles", choices=["neutron", "photon"],
                        nargs="+", default=["neutron", "photon"],
                        help="Incident particles to include")
    parser.add_argument("--temperatures", type=float,
                        default=[250.0, 293.6, 600.0, 900.0, 1200.0, 2500.0],
                        help="Temperatures in Kelvin", nargs="+")
    parser.add_argument("--nuclides", nargs="+", metavar="NUCLIDE",
                        help="Only convert these nuclides (e.g. Fe56 U235)")
    parser.add_argument("--cleanup", action="store_true",
                        help="Remove source files after conversion")
    parser.add_argument("--force", action="store_true",
                        help="Reconvert nuclides whose output already exists "
                             "(default: skip them to save time)")
    args = parser.parse_args()

    info = ENDF_RELEASES[args.release]
    endf_dir = find_or_download_endf(args.release, args.particles)

    if args.destination is None:
        args.destination = Path.home() / "nuclear_data" / info["dest"]

    lib_name = info["library"]
    print(f"Output: {args.destination}")
    print(f"Temperatures: {args.temperatures}")

    for particle in args.particles:
        (args.destination / particle).mkdir(parents=True, exist_ok=True)

    # Neutrons — parallel via NJOY
    if "neutron" in args.particles:
        endf_files = sorted(endf_dir.rglob("n-*.endf"))
        # Skip free neutron (no bound cross sections)
        endf_files = [f for f in endf_files if f.name != "n-000_n_001.endf"]
        endf_files = nuclide_filter(endf_files, args.nuclides)

        neutron_out = args.destination / "neutron"
        skipped = []
        if not args.force:
            todo = []
            for f in endf_files:
                name = _nuclide_name_from_neutron_endf(f)
                if name and (neutron_out / f"{name}.arrow" / "version.json").is_file():
                    skipped.append(name)
                else:
                    todo.append(f)
            endf_files = todo

        print(f"Found {len(endf_files) + len(skipped)} neutron ENDF files "
              f"({len(skipped)} already converted, skipping; use --force to reconvert)")

        failed = []
        total = len(endf_files)
        with Pool() as pool:
            results = []
            for f in endf_files:
                r = pool.apply_async(
                    convert_neutron,
                    (f, neutron_out),
                    dict(source_format="endf", temperatures=args.temperatures,
                         library=lib_name),
                )
                results.append((f, r))
            for i, (f, r) in enumerate(results, 1):
                try:
                    r.get()
                    print(f"[{i}/{total}] {f.stem}")
                except Exception as e:
                    print(f"[{i}/{total}] FAILED: {f.name}: {e}")
                    failed.append(f.name)

        if failed:
            print(f"\n{len(failed)} neutron files failed: {failed}")

    # Photons — sequential (fast, no NJOY)
    if "photon" in args.particles:
        photo_files, atom_files = find_photon_files(endf_dir)
        if args.nuclides:
            # For photon files, filter both lists in lockstep
            indices = [i for i, f in enumerate(photo_files)
                       if f in nuclide_filter([f], args.nuclides)]
            photo_files = [photo_files[i] for i in indices]
            atom_files = [atom_files[i] for i in indices]
        print(f"Found {len(photo_files)} photoatomic + {len(atom_files)} atomic relaxation files")
        total_photon = len(photo_files)
        for i, (photo_path, atom_path) in enumerate(zip(photo_files, atom_files), 1):
            print(f"[{i}/{total_photon}] {photo_path.stem}")
            convert_photon(
                photo_path, args.destination / "photon",
                atom_path=atom_path, library=lib_name,
            )

    write_index(args.destination)

    if args.cleanup:
        from shutil import rmtree
        print(f"Cleaning up source: {endf_dir}")
        rmtree(endf_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
