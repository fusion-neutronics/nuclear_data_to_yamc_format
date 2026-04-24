"""Download and extraction utilities for nuclear data files."""

import shutil
import ssl
import subprocess
from pathlib import Path
from urllib.request import Request, urlopen


def download_file(url, dest_dir, *, verify_ssl=True):
    """Download a file from *url* into *dest_dir*, skipping if already present.

    Writes to a ``.part`` sidecar first and atomically renames on success, so
    an interrupted download never leaves a truncated file that would fool a
    later run into thinking it had a complete cache.

    Returns the path to the downloaded file.
    """
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    filename = url.rsplit("/", 1)[-1]
    dest = dest_dir / filename
    if dest.exists():
        print(f"  Already downloaded: {dest}")
        return dest

    part = dest.with_suffix(dest.suffix + ".part")
    if part.exists():
        part.unlink()  # discard any stale partial from a prior crash

    print(f"  Downloading {url}")
    kwargs = {}
    if not verify_ssl:
        kwargs["context"] = ssl._create_unverified_context()
    req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    try:
        with urlopen(req, **kwargs) as resp, open(part, "wb") as f:
            total = int(resp.headers.get("Content-Length", 0))
            downloaded = 0
            while True:
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)
                if total:
                    pct = 100 * downloaded / total
                    print(
                        f"\r  {downloaded / 1e6:.0f}/{total / 1e6:.0f} MB ({pct:.0f}%)",
                        end="",
                        flush=True,
                    )
            if total:
                print()
        if total and downloaded != total:
            raise IOError(
                f"Download incomplete: got {downloaded} bytes, expected {total} "
                f"({url})"
            )
        part.replace(dest)
    except BaseException:
        if part.exists():
            part.unlink()
        raise
    return dest


def extract_archive(archive, dest_dir):
    """Extract a .zip, .tar.gz, .tar.xz, or single .endf file into *dest_dir*."""
    archive = Path(archive)
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    name = archive.name

    if name.endswith(".zip"):
        subprocess.run(
            ["unzip", "-o", str(archive), "-d", str(dest_dir)], check=True
        )
    elif name.endswith(".tar.gz") or name.endswith(".tgz"):
        subprocess.run(
            ["tar", "-xzf", str(archive), "-C", str(dest_dir)], check=True
        )
    elif name.endswith(".tar.xz"):
        subprocess.run(
            ["tar", "-xJf", str(archive), "-C", str(dest_dir)], check=True
        )
    elif name.endswith(".endf"):
        # Single ENDF file — copy into neutron subdir
        neutron_dir = dest_dir / "neutron"
        neutron_dir.mkdir(exist_ok=True)
        shutil.copy2(archive, neutron_dir / name)
    else:
        print(f"  Warning: unknown archive format: {name}")


def download_and_extract(urls, dest_dir, download_dir, *, verify_ssl=True):
    """Download a list of URLs and extract each archive into *dest_dir*.

    Parameters
    ----------
    urls : list of str
        Full URLs to download.
    dest_dir : Path
        Directory to extract files into.
    download_dir : Path
        Directory to store downloaded archives.
    verify_ssl : bool
        Whether to verify SSL certificates (some IAEA servers need False).
    """
    for url in urls:
        archive = download_file(url, download_dir, verify_ssl=verify_ssl)
        extract_archive(archive, dest_dir)


def find_photon_files(endf_dir):
    """Find photoatomic and atomic relaxation ENDF files anywhere under *endf_dir*.

    Matches files named ``photoat-*.endf`` and ``atom-*.endf`` regardless of
    the intermediate directory layout (e.g. ``photoat/``, ``photon/``,
    ``photoat-version.VIII.1/``).

    Returns ``(photo_files, atom_files)`` as sorted lists of Paths.
    """
    endf_dir = Path(endf_dir)
    return (
        sorted(endf_dir.rglob("photoat-*.endf")),
        sorted(endf_dir.rglob("atom-*.endf")),
    )


# =============================================================================
# Release metadata
# =============================================================================

ENDF_RELEASES = {
    "viii.1": {
        "library": "endfb-8.1",
        "dir": "endfb-viii.1-endf",
        "dest": "endf-b8.1-arrow",
        "neutron": {
            "base_url": "https://www.nndc.bnl.gov/endf-releases/releases/B-VIII.1/",
            "files": [
                "neutrons/neutrons-version.VIII.1.tar.gz",
            ],
        },
        "photon": {
            "base_url": "https://www.nndc.bnl.gov/endf-releases/releases/B-VIII.1/",
            "files": [
                "photoat/photoat-version.VIII.1.tar.gz",
                "atomic_relax/atomic_relax-version.VIII.1.tar.gz",
            ],
        },
        "decay": {
            "base_url": "https://www.nndc.bnl.gov/endf-releases/releases/B-VIII.1/",
            "files": [
                "decay/decay-version.VIII.1.tar.gz",
            ],
        },
        "nfy": {
            "base_url": "https://www.nndc.bnl.gov/endf-releases/releases/B-VIII.1/",
            "files": [
                "nfy/nfy-version.VIII.1.tar.gz",
            ],
        },
    },
    "viii.0": {
        "library": "endfb-8.0",
        "dir": "endfb-viii.0-endf",
        "dest": "endf-b8.0-arrow",
        "neutron": {
            "base_url": "https://www.nndc.bnl.gov/endf-b8.0/",
            "files": [
                "zips/ENDF-B-VIII.0_neutrons.zip",
                "zips/ENDF-B-VIII.0_thermal_scatt.zip",
                "erratafiles/n-005_B_010.endf",
            ],
        },
        "photon": {
            "base_url": "https://www.nndc.bnl.gov/endf-b8.0/",
            "files": [
                "zips/ENDF-B-VIII.0_photoat.zip",
                "erratafiles/atomic_relax.tar.gz",
            ],
        },
    },
    "vii.1": {
        "library": "endfb-7.1",
        "dir": "endfb-vii.1-endf",
        "dest": "endf-b7.1-arrow",
        "neutron": {
            "base_url": "http://www.nndc.bnl.gov/endf-b7.1/",
            "files": [
                "zips/ENDF-B-VII.1-neutrons.zip",
            ],
        },
        "photon": {
            "base_url": "http://www.nndc.bnl.gov/endf-b7.1/zips/",
            "files": [
                "ENDF-B-VII.1-photoat.zip",
                "ENDF-B-VII.1-atomic_relax.zip",
            ],
        },
    },
}

TENDL_RELEASES = {
    "2025": {
        "library": "tendl-2025",
        "dir": "tendl-2025-endf",
        "dest": "tendl-2025-arrow",
        "neutron": {
            "base_url": "https://tendl.imperial.ac.uk/tendl_2025/tar_files/",
            "files": ["TENDL-n.tgz"],
            "glob": "n-*.tendl",
        },
    },
    "2023": {
        "library": "tendl-2023",
        "dir": "tendl-2023-endf",
        "dest": "tendl-2023-arrow",
        "neutron": {
            "base_url": "https://tendl.imperial.ac.uk/tendl_2023/tar_files/",
            "files": ["TENDL-n.2024new.tgz"],
            "glob": "n-*.tendl",
        },
    },
}

FENDL_RELEASES = {
    "3.2c": {
        "library": "fendl-3.2c",
        "neutron": {
            "ace": {
                "base_url": "https://nds.iaea.org/fendl/data/neutron/",
                "files": ["fendl-FENDL-3.2c-neutron-ace.zip"],
                "glob": "neutron/ace/*",
            },
            "endf": {
                "base_url": "https://nds.iaea.org/fendl/data/neutron/",
                "files": ["fendl-FENDL-3.2c-neutron-endf.zip"],
                "glob": "neutron/endf/*.endf",
            },
        },
        "photon": {
            "endf": {
                "base_url": "https://nds.iaea.org/fendl/data/atom/",
                "files": ["fendl-FENDL-3.2c-atom-endf.zip"],
                "glob": "atom/endf/*.endf",
            },
        },
    },
    "3.1d": {
        "library": "fendl-3.1d",
        "neutron": {
            "ace": {
                "base_url": "https://nds.iaea.org/fendl/data/neutron/",
                "files": ["fendl-FENDL-3.1d-neutron-ace.zip"],
                "glob": "neutron/ace/*",
            },
        },
        "photon": {
            "endf": {
                "base_url": "https://nds.iaea.org/fendl/data/atom/",
                "files": ["fendl-FENDL-3.1d-atom-endf.zip"],
                "glob": "atom/endf/*.endf",
            },
        },
    },
}
