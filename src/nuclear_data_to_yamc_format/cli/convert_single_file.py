#!/usr/bin/env python3

"""Convert a single ACE or ENDF nuclear data file to simulation-ready Arrow format.

Examples:
    # Neutron from ACE
    convert-single-file neutron 92235.710nc

    # Neutron from ENDF (runs NJOY internally)
    convert-single-file neutron n-092_U_235.endf -f endf

    # Photon from ENDF
    convert-single-file photon photoat-026_Fe_000.endf

    # Photon from paired ENDF files
    convert-single-file photon photoat-026_Fe_000.endf --atom atom-026_Fe_000.endf

    # With library name
    convert-single-file neutron 92235.710nc --library endfb-8.0
"""

import argparse
from pathlib import Path

from nuclear_data_to_yamc_format import convert_neutron, convert_photon


def main():
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
    parser.add_argument("--library", type=str, default="",
                        help="Library name (e.g., endfb-8.0, fendl-3.2c)")
    args = parser.parse_args()

    output_dir = args.output or Path.cwd()

    if args.particle == "neutron":
        arrow_path = convert_neutron(
            args.input, output_dir,
            source_format=args.format,
            temperatures=args.temperatures,
            library=args.library,
        )
    else:
        arrow_path = convert_photon(
            args.input, output_dir,
            atom_path=args.atom,
            library=args.library,
        )

    print(f"Wrote {arrow_path}")


if __name__ == "__main__":
    main()
