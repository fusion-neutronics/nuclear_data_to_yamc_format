#!/usr/bin/env python3

"""Convert an OpenMC depletion chain to a simulation-ready .chain.arrow/ directory.

Examples:
    # Build from an existing OpenMC chain XML
    convert-chain --xml chain_endf_b8.0.xml \\
        -o chain_endf_b8.0.chain.arrow --library endfb-8.0

    # Build from ENDF source files (all reactions included, matching
    # openmc_data's generate_endf_chain.py)
    convert-chain \\
        --decay-dir decay/ --fpy-dir nfy/ --neutron-dir neutrons/ \\
        -o chain_endf_b8.0.chain.arrow --library endfb-8.0
"""

import argparse
from pathlib import Path

from nuclear_data_to_yamc_format import convert_chain


def _collect(dir_path):
    return sorted(Path(dir_path).glob("*.endf"))


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--xml", type=Path, default=None,
                        help="Existing OpenMC chain XML to convert")
    parser.add_argument("--decay-dir", type=Path, default=None,
                        help="Directory of decay ENDF files")
    parser.add_argument("--fpy-dir", type=Path, default=None,
                        help="Directory of neutron fission product yield ENDF files")
    parser.add_argument("--neutron-dir", type=Path, default=None,
                        help="Directory of neutron ENDF files")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="Output .chain.arrow/ directory")
    parser.add_argument("--library", type=str, default="",
                        help="Library name (e.g., endfb-8.0)")
    args = parser.parse_args()

    if args.xml is not None:
        path = convert_chain(args.output, xml_path=args.xml, library=args.library)
    else:
        if not all([args.decay_dir, args.fpy_dir, args.neutron_dir]):
            parser.error(
                "Either --xml or all of --decay-dir/--fpy-dir/--neutron-dir must be provided."
            )
        path = convert_chain(
            args.output,
            decay_files=_collect(args.decay_dir),
            fpy_files=_collect(args.fpy_dir),
            neutron_files=_collect(args.neutron_dir),
            library=args.library,
        )

    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
