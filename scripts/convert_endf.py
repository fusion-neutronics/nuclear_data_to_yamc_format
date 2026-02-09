#!/usr/bin/env python3

"""
Download ENDF/B-VIII.0 or ENDF/B-VII.1 library and convert to simulation-ready Arrow format.

Mirrors openmc_data generate_endf.py but replaces HDF5 export with Arrow export.
The ENDF→NJOY→ACE→OpenMC pipeline stays the same; only the final export changes.
"""

import argparse
import sys
import warnings
from multiprocessing import Pool
from pathlib import Path

import openmc.data
from openmc_data import download, state_download_size

from nuclear_data_to_yamc_format import convert_neutron, convert_photon
from nuclear_data_to_yamc_format.neutron_writer import export_neutron_to_arrow

assert sys.version_info >= (3, 9), "Python 3.9+ is required"


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass


parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=CustomFormatter
)
parser.add_argument('-d', '--destination', type=Path,
                    help='Directory to create new library in')
parser.add_argument('--download', action='store_true',
                    help='Download files')
parser.add_argument('--no-download', dest='download', action='store_false',
                    help='Do not download files')
parser.add_argument('--extract', action='store_true',
                    help='Extract zip files')
parser.add_argument('--no-extract', dest='extract', action='store_false',
                    help='Do not extract zip files')
parser.add_argument('-r', '--release', choices=['vii.1', 'viii.0'],
                    default='viii.0', help="ENDF/B release version")
parser.add_argument('-p', '--particles', choices=['neutron', 'photon'],
                    nargs='+', default=['neutron', 'photon'],
                    help="Incident particles to include")
parser.add_argument('--cleanup', action='store_true',
                    help="Remove download directories after processing")
parser.add_argument('--no-cleanup', dest='cleanup', action='store_false',
                    help="Do not remove download directories")
parser.add_argument('--temperatures', type=float,
                    default=[250.0, 293.6, 600.0, 900.0, 1200.0, 2500.0],
                    help="Temperatures in Kelvin", nargs='+')
parser.set_defaults(download=True, extract=True, cleanup=False)
args = parser.parse_args()

# Map release to library name
_LIBRARY_NAMES = {
    'viii.0': 'endfb-8.0',
    'vii.1': 'endfb-7.1',
}


def process_neutron_arrow(endf_path, output_dir, temperatures=None, library=""):
    """Process ENDF neutron file via NJOY and export to Arrow."""
    print(f'Converting: {endf_path}')
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            data = openmc.data.IncidentNeutron.from_njoy(
                endf_path, temperatures=temperatures
            )
    except Exception as e:
        print(f'{endf_path}: {e}')
        raise
    arrow_dir = output_dir / f'{data.name}.yamc.arrow'
    print(f'Writing {arrow_dir} ...')
    export_neutron_to_arrow(data, arrow_dir, library=library)


def main():
    library_name = 'endfb'
    cwd = Path.cwd()

    endf_files_dir = cwd / '-'.join([library_name, args.release, 'endf'])
    neutron_dir = endf_files_dir / 'neutron'
    download_path = cwd / '-'.join([library_name, args.release, 'download'])

    if args.destination is None:
        args.destination = Path('-'.join([library_name, args.release, 'arrow']))

    lib_name = _LIBRARY_NAMES.get(args.release, f'endfb-{args.release}')

    # Release details (same URLs as openmc_data generate_endf.py)
    release_details = {
        'viii.0': {
            'neutron': {
                'base_url': 'https://www.nndc.bnl.gov/endf-b8.0/',
                'compressed_files': ['zips/ENDF-B-VIII.0_neutrons.zip',
                                     'zips/ENDF-B-VIII.0_thermal_scatt.zip',
                                     'erratafiles/n-005_B_010.endf'],
                'endf_files': neutron_dir.rglob('n-*.endf'),
            },
            'photon': {
                'base_url': 'https://www.nndc.bnl.gov/endf-b8.0/',
                'compressed_files': ['zips/ENDF-B-VIII.0_photoat.zip',
                                     'erratafiles/atomic_relax.tar.gz'],
                'photo_files': endf_files_dir.joinpath('photoat').rglob('*.endf'),
                'atom_files': endf_files_dir.joinpath('atom').rglob('*.endf'),
            },
        },
    }

    # Create output directories
    for particle in args.particles:
        particle_destination = args.destination / particle
        particle_destination.mkdir(parents=True, exist_ok=True)

    # =========================================================================
    # PROCESS INCIDENT NEUTRON DATA
    if 'neutron' in args.particles:
        particle = 'neutron'
        details = release_details[args.release][particle]
        with Pool() as pool:
            results = []
            for filename in details['endf_files']:
                if filename.name == 'n-000_n_001.endf':
                    continue
                func_args = (filename, args.destination / particle,
                             args.temperatures, lib_name)
                r = pool.apply_async(process_neutron_arrow, func_args)
                results.append(r)
            for r in results:
                r.wait()

    # =========================================================================
    # INCIDENT PHOTON DATA
    if 'photon' in args.particles:
        particle = 'photon'
        details = release_details[args.release][particle]
        for photo_path, atom_path in zip(sorted(details['photo_files']),
                                         sorted(details['atom_files'])):
            print('Converting:', photo_path.name, atom_path.name)
            convert_photon(photo_path, args.destination / particle,
                           atom_path=atom_path, library=lib_name)


if __name__ == '__main__':
    main()
