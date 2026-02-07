#!/usr/bin/env python3

"""Convert a single ACE or ENDF nuclear data file to Arrow format.

Examples:
    # Neutron from ACE
    python scripts/convert_single_file.py neutron 92235.710nc

    # Neutron from ENDF (runs NJOY internally)
    python scripts/convert_single_file.py neutron n-092_U_235.endf -f endf

    # Photon from ENDF
    python scripts/convert_single_file.py photon photoat-026_Fe_000.endf

    # Photon from paired ENDF files
    python scripts/convert_single_file.py photon photoat-026_Fe_000.endf --atom atom-026_Fe_000.endf
"""

import argparse
from pathlib import Path

from nuclear_data_to_yamc_format import convert_neutron, convert_photon

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("particle", choices=["neutron", "photon"])
parser.add_argument("input", type=Path, help="Path to an ACE or ENDF file")
parser.add_argument("--atom", type=Path, default=None,
                    help="Atomic relaxation ENDF file (photon only)")
parser.add_argument("-f", "--format", choices=["ace", "endf"], default="ace",
                    help="Input file format (default: ace)")
parser.add_argument("-o", "--output", type=Path, default=None,
                    help="Output directory (default: current directory)")
parser.add_argument("--temperatures", type=float, nargs="+", default=None,
                    help="Temperatures in Kelvin (ENDF neutron only)")
args = parser.parse_args()


def main():
    output_dir = args.output or Path.cwd()

    if args.particle == "neutron":
        arrow_path = convert_neutron(
            args.input, output_dir,
            source_format=args.format,
            temperatures=args.temperatures,
        )
    else:
        arrow_path = convert_photon(
            args.input, output_dir,
            atom_path=args.atom,
        )

    print(f"Wrote {arrow_path}")


if __name__ == "__main__":
    main()
