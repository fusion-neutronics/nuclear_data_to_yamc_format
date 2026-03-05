#!/usr/bin/env python3

"""
Download and convert ENDF/B-VIII.0 or VII.1 to simulation-ready Arrow format.

Pipeline: download ENDF → NJOY (Doppler broaden) → ACE → OpenMC parse → Arrow.

If ENDF source files already exist locally, the download step is skipped.
Search order: ./endfb-{release}-endf/ then ~/nuclear_data/endfb-{release}-endf/.
Output defaults to ~/nuclear_data/endf-b{X}.0-arrow/.
"""

import argparse
import sys
import warnings
from multiprocessing import Pool
from pathlib import Path

import openmc.data

from nuclear_data_to_yamc_format import convert_photon
from nuclear_data_to_yamc_format.download import (
    ENDF_RELEASES, download_and_extract, find_photon_files,
)
from nuclear_data_to_yamc_format.neutron_writer import export_neutron_to_arrow

assert sys.version_info >= (3, 9), "Python 3.9+ is required"


def process_neutron_arrow(endf_path, output_dir, temperatures=None, library=""):
    """Process ENDF neutron file via NJOY and export to Arrow."""
    print(f"Converting: {endf_path}")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            data = openmc.data.IncidentNeutron.from_njoy(
                endf_path, temperatures=temperatures
            )
    except Exception as e:
        print(f"{endf_path}: {e}")
        raise
    arrow_dir = output_dir / f"{data.name}.arrow"
    print(f"Writing {arrow_dir} ...")
    export_neutron_to_arrow(data, arrow_dir, library=library)


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
    args = parser.parse_args()

    info = ENDF_RELEASES[args.release]
    endf_dir = find_or_download_endf(args.release, args.particles)
    neutron_dir = endf_dir / "neutron"

    if args.destination is None:
        args.destination = Path.home() / "nuclear_data" / info["dest"]

    lib_name = info["library"]
    print(f"Output: {args.destination}")
    print(f"Temperatures: {args.temperatures}")

    for particle in args.particles:
        (args.destination / particle).mkdir(parents=True, exist_ok=True)

    # Neutrons
    if "neutron" in args.particles:
        endf_files = sorted(neutron_dir.rglob("n-*.endf"))
        print(f"Found {len(endf_files)} neutron ENDF files")
        with Pool() as pool:
            results = []
            for f in endf_files:
                if f.name == "n-000_n_001.endf":
                    continue
                r = pool.apply_async(
                    process_neutron_arrow,
                    (f, args.destination / "neutron", args.temperatures, lib_name),
                )
                results.append(r)
            for r in results:
                r.wait()

    # Photons
    if "photon" in args.particles:
        photo_files, atom_files = find_photon_files(endf_dir)
        print(f"Found {len(photo_files)} photoatomic + {len(atom_files)} atomic relaxation files")
        for photo_path, atom_path in zip(photo_files, atom_files):
            print("Converting:", photo_path.name, atom_path.name)
            convert_photon(
                photo_path, args.destination / "photon",
                atom_path=atom_path, library=lib_name,
            )


if __name__ == "__main__":
    main()
