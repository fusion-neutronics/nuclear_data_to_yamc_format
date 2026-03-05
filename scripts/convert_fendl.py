#!/usr/bin/env python3

"""
Download FENDL ACE/ENDF files and convert to simulation-ready Arrow format.
"""

import argparse
from pathlib import Path
from shutil import rmtree

import openmc.data

from nuclear_data_to_yamc_format import convert_neutron
from nuclear_data_to_yamc_format.download import (
    FENDL_RELEASES, download_and_extract,
)
from nuclear_data_to_yamc_format.photon_writer import export_photon_to_arrow

# Default file format per particle type
FILE_TYPES = {"neutron": "ace", "photon": "endf"}


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
    parser.add_argument("--download", action="store_true", help="Download files")
    parser.add_argument("--no-download", dest="download", action="store_false",
                        help="Do not download files")
    parser.add_argument("--extract", action="store_true", help="Extract tar/zip files")
    parser.add_argument("--no-extract", dest="extract", action="store_false",
                        help="Do not extract tar/zip files")
    parser.add_argument("-r", "--release", choices=list(FENDL_RELEASES.keys()),
                        default="3.2c", help="FENDL release version")
    parser.add_argument("-p", "--particles", choices=["neutron", "photon"],
                        nargs="+", default=["neutron", "photon"],
                        help="Incident particles to include")
    parser.add_argument("--cleanup", action="store_true",
                        help="Remove download directories after processing")
    parser.add_argument("--no-cleanup", dest="cleanup", action="store_false",
                        help="Do not remove download directories")
    parser.set_defaults(download=True, extract=True, cleanup=False)
    args = parser.parse_args()

    library_name = "fendl"
    cwd = Path.cwd()
    release = FENDL_RELEASES[args.release]
    lib_name = release["library"]

    ace_files_dir = cwd / f"{library_name}-{args.release}-ace"
    endf_files_dir = cwd / f"{library_name}-{args.release}-endf"
    download_dir = cwd / f"{library_name}-{args.release}-download"

    if args.destination is None:
        args.destination = Path(f"{library_name}-{args.release}-arrow")

    # ==========================================================================
    # DOWNLOAD & EXTRACT
    if args.download or args.extract:
        for particle in args.particles:
            file_type = FILE_TYPES[particle]
            details = release[particle][file_type]
            urls = [details["base_url"] + f for f in details["files"]]
            extraction_dir = ace_files_dir if file_type == "ace" else endf_files_dir

            if args.download:
                # IAEA servers sometimes have cert issues
                download_and_extract(
                    urls, extraction_dir, download_dir / particle,
                    verify_ssl=False,
                )
            elif args.extract:
                # Extract only (already downloaded)
                from nuclear_data_to_yamc_format.download import extract_archive
                for url in urls:
                    filename = url.rsplit("/", 1)[-1]
                    archive = download_dir / particle / filename
                    if archive.exists():
                        extract_archive(archive, extraction_dir)

        if args.cleanup and download_dir.exists():
            rmtree(download_dir)

    # ==========================================================================
    # GENERATE ARROW LIBRARY
    for particle in args.particles:
        particle_dest = args.destination / particle
        particle_dest.mkdir(parents=True, exist_ok=True)

        file_type = FILE_TYPES[particle]
        details = release[particle][file_type]

        if particle == "neutron":
            neutron_files = ace_files_dir.glob(details["glob"])
            neutron_files = [
                f for f in neutron_files
                if not f.name.endswith("_") and not f.name.endswith(".xsd")
            ]
            for filename in sorted(neutron_files):
                print(f"Converting: {filename}")
                convert_neutron(filename, particle_dest,
                                source_format="ace", library=lib_name)

            if args.cleanup and ace_files_dir.exists():
                rmtree(ace_files_dir)

        elif particle == "photon":
            photon_files = endf_files_dir.glob(details["glob"])
            for photo_path in sorted(photon_files):
                print(f"Converting: {photo_path}")
                evaluations = openmc.data.endf.get_evaluations(photo_path)
                for ev in evaluations:
                    data = openmc.data.IncidentPhoton.from_endf(ev)
                    arrow_dir = particle_dest / f"{data.name}.arrow"
                    print(f"Writing {arrow_dir}...")
                    export_photon_to_arrow(data, arrow_dir, library=lib_name)

            if args.cleanup and endf_files_dir.exists():
                rmtree(endf_files_dir)


if __name__ == "__main__":
    main()
