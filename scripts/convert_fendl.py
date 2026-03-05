#!/usr/bin/env python3

"""
Convert FENDL to simulation-ready Arrow format.

Downloads ACE (neutron) and ENDF (photon) files from IAEA if not found locally.
Search order: ./fendl-{release}-{ace,endf}/ then ~/nuclear_data/fendl-{release}-{ace,endf}/.
Output defaults to ./fendl-{release}-arrow/.
"""

import argparse
from pathlib import Path

from nuclear_data_to_yamc_format import convert_neutron, convert_photon_endf
from nuclear_data_to_yamc_format.download import (
    FENDL_RELEASES, download_and_extract, extract_archive,
)

# Default file format per particle type
FILE_TYPES = {"neutron": "ace", "photon": "endf"}


def find_or_download_fendl(release_key, particles):
    """Find existing FENDL source files or download them.

    Returns (ace_dir, endf_dir) — either may be None if that particle
    type was not requested.
    """
    release = FENDL_RELEASES[release_key]
    lib_name = f"fendl-{release_key}"

    ace_dir = None
    endf_dir = None

    for particle in particles:
        file_type = FILE_TYPES[particle]
        details = release[particle][file_type]
        suffix = "ace" if file_type == "ace" else "endf"
        dirname = f"{lib_name}-{suffix}"

        # Search for existing files
        found = None
        for base in [Path.cwd(), Path.home() / "nuclear_data"]:
            candidate = base / dirname
            if candidate.is_dir() and any(candidate.rglob("*")):
                found = candidate
                break

        if found is None:
            # Download
            found = Path.cwd() / dirname
            download_dir = Path.cwd() / f"{lib_name}-download" / particle
            urls = [details["base_url"] + f for f in details["files"]]
            print(f"Downloading {particle} ({file_type}) data...")
            download_and_extract(urls, found, download_dir, verify_ssl=False)

        print(f"Using {particle} source: {found}")
        if file_type == "ace":
            ace_dir = found
        else:
            endf_dir = found

    return ace_dir, endf_dir


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
    parser.add_argument("-r", "--release", choices=list(FENDL_RELEASES.keys()),
                        default="3.2c", help="FENDL release version")
    parser.add_argument("-p", "--particles", choices=["neutron", "photon"],
                        nargs="+", default=["neutron", "photon"],
                        help="Incident particles to include")
    parser.add_argument("--cleanup", action="store_true",
                        help="Remove source files after conversion")
    args = parser.parse_args()

    release = FENDL_RELEASES[args.release]
    lib_name = release["library"]
    ace_dir, endf_dir = find_or_download_fendl(args.release, args.particles)

    if args.destination is None:
        args.destination = Path(f"fendl-{args.release}-arrow")

    # Neutrons — from ACE
    if "neutron" in args.particles:
        neutron_dest = args.destination / "neutron"
        neutron_dest.mkdir(parents=True, exist_ok=True)

        details = release["neutron"]["ace"]
        neutron_files = sorted(ace_dir.glob(details["glob"]))
        neutron_files = [
            f for f in neutron_files
            if not f.name.endswith("_") and not f.name.endswith(".xsd")
        ]
        print(f"Found {len(neutron_files)} neutron ACE files")

        for filename in neutron_files:
            print(f"Converting: {filename.name}")
            convert_neutron(filename, neutron_dest,
                            source_format="ace", library=lib_name)

    # Photons — from ENDF (may have multiple evaluations per file)
    if "photon" in args.particles:
        photon_dest = args.destination / "photon"
        photon_dest.mkdir(parents=True, exist_ok=True)

        details = release["photon"]["endf"]
        photon_files = sorted(endf_dir.glob(details["glob"]))
        print(f"Found {len(photon_files)} photon ENDF files")

        for photo_path in photon_files:
            print(f"Converting: {photo_path.name}")
            convert_photon_endf(photo_path, photon_dest, library=lib_name)

    if args.cleanup:
        from shutil import rmtree
        if ace_dir:
            print(f"Cleaning up: {ace_dir}")
            rmtree(ace_dir, ignore_errors=True)
        if endf_dir:
            print(f"Cleaning up: {endf_dir}")
            rmtree(endf_dir, ignore_errors=True)
        dl_dir = Path.cwd() / f"fendl-{args.release}-download"
        if dl_dir.exists():
            print(f"Cleaning up: {dl_dir}")
            rmtree(dl_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
