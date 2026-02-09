from importlib.metadata import version, PackageNotFoundError
from pathlib import Path

import openmc.data

from .neutron_writer import export_neutron_to_arrow
from .photon_writer import export_photon_to_arrow
from .neutron_reader import read_neutron_from_arrow
from .photon_reader import read_photon_from_arrow
from .verify import verify_neutron, verify_photon

try:
    __version__ = version("nuclear_data_to_yamc_format")
except PackageNotFoundError:
    __version__ = "0.0.0-dev"


def convert_neutron(input_path, output_dir, *, source_format="ace",
                    temperatures=None, library=""):
    """Convert an ACE or ENDF neutron file to simulation-ready Arrow format.

    Parameters
    ----------
    input_path : str or Path
        Path to an ACE or ENDF neutron file.
    output_dir : str or Path
        Directory to write the Arrow output into.
        A subdirectory named ``{nuclide}.arrow/`` will be created.
    source_format : {"ace", "endf"}
        Input file format. ``"endf"`` runs NJOY internally.
    temperatures : list of float, optional
        Temperatures in Kelvin (only used with ``source_format="endf"``).
    library : str, optional
        Library name (e.g., "endfb-8.0", "fendl-3.2c").

    Returns
    -------
    Path
        Path to the created Arrow directory.
    """
    input_path = Path(input_path)
    output_dir = Path(output_dir)

    if source_format == "ace":
        data = openmc.data.IncidentNeutron.from_ace(input_path)
    elif source_format == "endf":
        kwargs = {}
        if temperatures is not None:
            kwargs["temperatures"] = temperatures
        data = openmc.data.IncidentNeutron.from_njoy(input_path, **kwargs)
    else:
        raise ValueError(f"Unknown source_format: {source_format!r}")

    arrow_path = output_dir / f"{data.name}.arrow"
    export_neutron_to_arrow(data, arrow_path, library=library)
    return arrow_path


def convert_photon(input_path, output_dir, *, atom_path=None, library=""):
    """Convert an ENDF photon file to simulation-ready Arrow format.

    Parameters
    ----------
    input_path : str or Path
        Path to a photoatomic ENDF file.
    output_dir : str or Path
        Directory to write the Arrow output into.
        A subdirectory named ``{element}.arrow/`` will be created.
    atom_path : str or Path, optional
        Path to an atomic relaxation ENDF file.
    library : str, optional
        Library name (e.g., "endfb-8.0", "fendl-3.2c").

    Returns
    -------
    Path
        Path to the created Arrow directory.
    """
    input_path = Path(input_path)
    output_dir = Path(output_dir)

    if atom_path is not None:
        data = openmc.data.IncidentPhoton.from_endf(input_path, atom_path)
    else:
        data = openmc.data.IncidentPhoton.from_endf(input_path)

    arrow_path = output_dir / f"{data.name}.arrow"
    export_photon_to_arrow(data, arrow_path, library=library)
    return arrow_path


__all__ = [
    "__version__",
    "convert_neutron",
    "convert_photon",
    "export_neutron_to_arrow",
    "export_photon_to_arrow",
    "read_neutron_from_arrow",
    "read_photon_from_arrow",
    "verify_neutron",
    "verify_photon",
]
