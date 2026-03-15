#!/usr/bin/env python3

"""
Convert TENDL to simulation-ready Arrow format.

Pipeline: download ENDF -> NJOY (Doppler broaden) -> OpenMC parse -> Arrow.

If ENDF source files already exist locally, the download step is skipped.
Search order: ./tendl-{release}-endf/ then ~/nuclear_data/tendl-{release}-endf/.
Output defaults to ./tendl-{release}-arrow/.

Only neutron data is available in TENDL.
"""

import argparse
import sys
from multiprocessing import Pool
from pathlib import Path

from nuclear_data_to_yamc_format import convert_neutron
from nuclear_data_to_yamc_format.cli import nuclide_filter
from nuclear_data_to_yamc_format.download import (
    TENDL_RELEASES, download_and_extract,
)

assert sys.version_info >= (3, 9), "Python 3.9+ is required"


def find_or_download_tendl(release, info):
    """Find existing TENDL ENDF source or download it."""
    dirname = info["dir"]

    candidates = [
        Path.cwd() / dirname,
        Path.home() / "nuclear_data" / dirname,
    ]
    for p in candidates:
        if p.is_dir() and any(p.glob("n-*.tendl")):
            print(f"Using existing TENDL source: {p}")
            return p

    dest = Path.cwd() / dirname
    download_dir = Path.cwd() / f"tendl-{release}-download"
    print(f"TENDL source not found locally. Downloading to {dest}")

    neutron = info["neutron"]
    urls = [neutron["base_url"] + f for f in neutron["files"]]
    print("Downloading neutron data...")
    download_and_extract(urls, dest, download_dir, verify_ssl=False)

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
    parser.add_argument("-d", "--destination", type=Path, default=None,
                        help="Directory to create new library in")
    parser.add_argument("-r", "--release", choices=list(TENDL_RELEASES.keys()),
                        default="2025", help="TENDL release version")
    parser.add_argument("--temperatures", type=float,
                        default=[250.0, 293.6, 600.0, 900.0, 1200.0, 2500.0],
                        help="Temperatures in Kelvin", nargs="+")
    parser.add_argument("--nuclides", nargs="+", metavar="NUCLIDE",
                        help="Only convert these nuclides (e.g. Fe56 U235)")
    parser.add_argument("--cleanup", action="store_true",
                        help="Remove source files after conversion")
    args = parser.parse_args()

    info = TENDL_RELEASES[args.release]
    lib_name = info["library"]
    endf_dir = find_or_download_tendl(args.release, info)

    if args.destination is None:
        args.destination = Path(info["dest"])

    neutron_dest = args.destination / "neutron"
    neutron_dest.mkdir(parents=True, exist_ok=True)

    print(f"Output: {args.destination}")
    print(f"Temperatures: {args.temperatures}")

    neutron_glob = info["neutron"]["glob"]
    endf_files = sorted(endf_dir.glob(neutron_glob))
    endf_files = nuclide_filter(endf_files, args.nuclides)
    print(f"Found {len(endf_files)} neutron ENDF files")

    failed = []
    total = len(endf_files)
    with Pool() as pool:
        results = []
        for f in endf_files:
            r = pool.apply_async(
                convert_neutron,
                (f, neutron_dest),
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
        print(f"\n{len(failed)} files failed: {failed}")

    if args.cleanup:
        from shutil import rmtree
        print(f"Cleaning up source: {endf_dir}")
        rmtree(endf_dir, ignore_errors=True)
        dl_dir = Path.cwd() / f"tendl-{args.release}-download"
        if dl_dir.exists():
            print(f"Cleaning up: {dl_dir}")
            rmtree(dl_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
