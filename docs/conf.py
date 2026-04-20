import importlib.metadata

project = "nuclear_data_to_yamc_format"
copyright = "2026, nuclear_data_to_yamc_format contributors"
author = "nuclear_data_to_yamc_format contributors"

try:
    release = importlib.metadata.version("nuclear_data_to_yamc_format")
except importlib.metadata.PackageNotFoundError:
    release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "myst_parser",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

myst_enable_extensions = [
    "colon_fence",
    "fieldlist",
    "deflist",
]

templates_path = ["_templates"]
exclude_patterns = ["_build"]

html_theme = "pydata_sphinx_theme"

html_theme_options = {
    "github_url": "https://github.com/yamc-org/nuclear_data_to_yamc_format",
    "show_toc_level": 2,
    "navigation_with_keys": False,
    "navbar_align": "left",
    "secondary_sidebar_items": ["page-toc", "edit-this-page"],
    "icon_links": [],
}

html_context = {
    "default_mode": "light",
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pyarrow": ("https://arrow.apache.org/docs/", None),
}

autodoc_member_order = "bysource"
autodoc_typehints = "description"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
