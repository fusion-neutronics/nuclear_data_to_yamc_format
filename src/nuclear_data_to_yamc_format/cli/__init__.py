import re
import sys
from pathlib import Path


def write_index(library_dir):
    """Write a comma-delimited index of all converted nuclides and elements.

    Scans ``neutron/`` and ``photon/`` subdirectories for ``*.arrow``
    directories and writes their names to ``index.txt`` in *library_dir*.
    """
    names = set()
    library_dir = Path(library_dir)

    for subdir in ("neutron", "photon"):
        particle_dir = library_dir / subdir
        if not particle_dir.is_dir():
            continue
        for d in particle_dir.iterdir():
            if d.is_dir() and d.suffix == ".arrow" and ".photon" not in d.stem:
                names.add(d.stem)

    if names:
        index_path = library_dir / "index.txt"
        index_path.write_text(",".join(sorted(names)) + "\n")
        print(f"Wrote index ({len(names)} entries) to {index_path}")


def parse_nuclide(name):
    """Parse 'Fe56' -> ('Fe', 56). Expects format like Fe56, U235, H1."""
    m = re.match(r"^([A-Z][a-z]?)(\d+)$", name)
    if m:
        return m.group(1), int(m.group(2))
    print(f"Error: invalid nuclide format {name!r}, expected e.g. Fe56 or U235")
    sys.exit(1)


def nuclide_filter(files, nuclides):
    """Keep only files whose nuclide is in *nuclides* (e.g. ['Fe56', 'U235']).

    Supports ENDF-style names (n-026_Fe_056.endf) and ACE-style names.
    Returns all files if *nuclides* is None or empty.
    Exits with an error if a requested nuclide has no matching file.
    """
    if not nuclides:
        return files

    targets = {parse_nuclide(n) for n in nuclides}

    kept = []
    matched = set()
    for f in files:
        stem = f.stem
        # ENDF / TENDL pattern: n-ZZZ_Sy_AAA  or  photoat-ZZZ_Sy_AAA
        parts = stem.split("-", 1)[-1].split("_")
        if len(parts) >= 3:
            try:
                elem, mass = parts[1], int(parts[2])
                if (elem, mass) in targets:
                    kept.append(f)
                    matched.add((elem, mass))
                continue
            except (ValueError, IndexError):
                pass

        # Fallback: check if any target nuclide appears in the filename
        # Try both unpadded (Fe56) and zero-padded (Fe056) mass numbers
        for elem, mass in targets:
            if f"{elem}{mass}" in stem or f"{elem}{mass:03d}" in stem:
                kept.append(f)
                matched.add((elem, mass))
                break

    missing = targets - matched
    if missing:
        names = ", ".join(f"{e}{m}" for e, m in sorted(missing))
        print(f"Error: no matching files found for: {names}")
        sys.exit(1)

    return kept
