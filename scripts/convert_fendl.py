#!/usr/bin/env python3

"""
Download FENDL ACE/ENDF files and convert to simulation-ready Arrow format.

Mirrors openmc_data convert_fendl.py but replaces HDF5 export with Arrow export.
"""

import argparse
import ssl
import subprocess
import warnings
from pathlib import Path
from shutil import rmtree
from urllib.parse import urljoin

import openmc.data
from openmc_data import download, all_release_details, calculate_download_size, get_file_types

from nuclear_data_to_yamc_format import convert_neutron, convert_photon
from nuclear_data_to_yamc_format.photon_writer import export_photon_to_arrow


class CustomFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter
):
    pass


parser = argparse.ArgumentParser(description=__doc__, formatter_class=CustomFormatter)
parser.add_argument("-d", "--destination", type=Path, default=None,
                    help="Directory to create new library in")
parser.add_argument("--download", action="store_true", help="Download files")
parser.add_argument("--no-download", dest="download", action="store_false",
                    help="Do not download files")
parser.add_argument("--extract", action="store_true", help="Extract tar/zip files")
parser.add_argument("--no-extract", dest="extract", action="store_false",
                    help="Do not extract tar/zip files")
parser.add_argument("-r", "--release",
                    choices=["3.2c", "3.2b", "3.2a", "3.2", "3.1d", "3.1a",
                             "3.1", "3.0", "2.1"],
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


def main():
    library_name = "fendl"
    file_types = get_file_types(args.particles)
    cwd = Path.cwd()

    ace_files_dir = cwd / "-".join([library_name, args.release, "ace"])
    endf_files_dir = cwd / "-".join([library_name, args.release, "endf"])

    download_path = cwd / "-".join([library_name, args.release, "download"])

    if args.destination is None:
        args.destination = Path("-".join([library_name, args.release, "arrow"]))

    lib_name = f"fendl-{args.release}"

    release_details = all_release_details[library_name]

    output_warnings = []

    # ==============================================================================
    # DOWNLOAD
    if args.download:
        calculate_download_size(library_name, args.release, args.particles,
                                file_types, 'GB')
        for particle in args.particles:
            particle_details = release_details[args.release][particle][file_types[particle]]
            for f in particle_details["compressed_files"]:
                download(
                    urljoin(particle_details["base_url"], f),
                    as_browser=True,
                    context=ssl._create_unverified_context(),
                    output_path=download_path / particle,
                )

    # ==============================================================================
    # EXTRACT
    if args.extract:
        for particle in args.particles:
            particle_details = release_details[args.release][particle][file_types[particle]]
            if file_types[particle] == "ace":
                extraction_dir = ace_files_dir
            elif file_types[particle] == "endf":
                extraction_dir = endf_files_dir

            for f in particle_details["compressed_files"]:
                subprocess.call(
                    ["unzip", "-o", download_path / particle / f, "-d",
                     extraction_dir]
                )

        if args.cleanup and download_path.exists():
            rmtree(download_path)

    # ==============================================================================
    # GENERATE ARROW LIBRARY

    for particle in args.particles:
        particle_destination = args.destination / particle
        particle_destination.mkdir(parents=True, exist_ok=True)

        particle_details = release_details[args.release][particle][file_types[particle]]

        if particle == "neutron":
            neutron_files = ace_files_dir.glob(
                release_details[args.release]["neutron"][file_types[particle]]["ace_files"]
            )
            neutron_files = [
                f for f in neutron_files
                if not f.name.endswith("_") and not f.name.endswith(".xsd")
            ]

            for filename in sorted(neutron_files):
                print(f"Converting: {filename}")
                convert_neutron(filename, particle_destination,
                                source_format="ace", library=lib_name)

            if args.cleanup and ace_files_dir.exists():
                rmtree(ace_files_dir)

        elif particle == "photon":
            photon_files = endf_files_dir.glob(
                release_details[args.release]["photon"][file_types[particle]]["endf_files"]
            )

            for photo_path in sorted(photon_files):
                print(f"Converting: {photo_path}")
                evaluations = openmc.data.endf.get_evaluations(photo_path)
                for ev in evaluations:
                    data = openmc.data.IncidentPhoton.from_endf(ev)
                    arrow_dir = particle_destination / f"{data.name}.photon.yamc.arrow"
                    print(f"Writing {arrow_dir}...")
                    export_photon_to_arrow(data, arrow_dir, library=lib_name)

            if args.cleanup and endf_files_dir.exists():
                rmtree(endf_files_dir)

    for warning in output_warnings:
        warnings.warn(warning)


if __name__ == "__main__":
    main()
