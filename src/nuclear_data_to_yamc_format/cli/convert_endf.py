#!/usr/bin/env python3

"""
Convert ENDF/B-VIII.0 or VII.1 to simulation-ready Arrow format.

Pipeline: download ENDF -> NJOY (Doppler broaden) -> OpenMC parse -> Arrow.

If ENDF source files already exist locally, the download step is skipped.
Search order: ./endfb-{release}-endf/ then ~/nuclear_data/endfb-{release}-endf/.
Output defaults to ~/nuclear_data/endf-b{X}.0-arrow/.
"""

import argparse
import sys
from multiprocessing import Pool
from pathlib import Path

from nuclear_data_to_yamc_format import convert_neutron, convert_photon
from nuclear_data_to_yamc_format.cli import nuclide_filter, write_index
from nuclear_data_to_yamc_format.download import (
    ENDF_RELEASES, download_and_extract, find_photon_files,
)

assert sys.version_info >= (3, 9), "Python 3.9+ is required"


def find_or_download_endf(release, particles):
    """Find existing ENDF source or download it."""
    info = ENDF_RELEASES[release]
    dirname = info["dir"]

    candidates = [
        Path.cwd() / dirname,
        Path.home() / "nuclear_data" / dirname,
    ]
    for p in candidates:
        if p.is_dir() and any(p.rglob("*.endf")):
            print(f"Using existing ENDF source: {p}")
            return p

    dest = Path.home() / "nuclear_data" / dirname
    download_dir = dest / "_downloads"
    print(f"ENDF source not found locally. Downloading to {dest}")

    for particle in particles:
        if particle in info:
            details = info[particle]
            urls = [details["base_url"] + f for f in details["files"]]
            print(f"\nDownloading {particle} data...")
            download_and_extract(urls, dest, download_dir)

    return dest


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
                        default="viii.0", help="ENDF/B release version")
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
        neutron_dir = endf_dir / "neutron"
        endf_files = sorted(neutron_dir.rglob("n-*.endf"))
        # Skip free neutron (no bound cross sections)
        endf_files = [f for f in endf_files if f.name != "n-000_n_001.endf"]
        endf_files = nuclide_filter(endf_files, args.nuclides)
        print(f"Found {len(endf_files)} neutron ENDF files")

        failed = []
        total = len(endf_files)
        with Pool() as pool:
            results = []
            for f in endf_files:
                r = pool.apply_async(
                    convert_neutron,
                    (f, args.destination / "neutron"),
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
