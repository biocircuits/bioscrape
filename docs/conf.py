from __future__ import annotations

import sys
import sysconfig
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Prefer an existing local build directory only when it contains extension
# modules compatible with the Python interpreter running Sphinx. This avoids
# accidentally importing stale cpXY artifacts from another local build. When
# RTD has installed the package into site-packages, keep the source checkout
# from shadowing the installed extension modules.
extension_suffix = sysconfig.get_config_var("EXT_SUFFIX")

for entry in list(sys.path):
    entry_path = Path.cwd() if entry == "" else Path(entry)
    try:
        if entry_path.resolve() == ROOT:
            sys.path.remove(entry)
    except OSError:
        pass

source_package = ROOT / "bioscrape"
if extension_suffix and any(source_package.glob(f"*{extension_suffix}")):
    sys.path.insert(0, str(ROOT))
else:
    sys.path.append(str(ROOT))

for build_lib in sorted((ROOT / "build").glob("lib*"), reverse=True):
    package_dir = build_lib / "bioscrape"
    if extension_suffix and any(package_dir.glob(f"*{extension_suffix}")):
        sys.path.insert(0, str(build_lib))

project = "bioscrape"
author = "Biocircuits, California Institute of Technology"
copyright = "2026, Biocircuits, California Institute of Technology"
release = "1.3.0"
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "nbsphinx",
]

autosummary_generate = True
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

autodoc_mock_imports = [
    "roadrunner",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "**.ipynb_checkpoints", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_title = "bioscrape documentation"

nbsphinx_execute = "never"
nbsphinx_allow_errors = True
nbsphinx_kernel_name = "python3"

# The examples are tutorial artifacts; their output cells are preserved and
# source links are more useful than long execution traces in RTD builds.
nbsphinx_timeout = 120


def _normalize_legacy_docstrings(app, what, name, obj, options, lines):
    """Make older Cython docstrings acceptable to docutils.

    Several public Cython methods use indented ``field: description`` blocks
    without the blank lines reStructuredText expects. Converting those lines to
    bullets keeps the content visible while avoiding noisy RTD warnings.
    """
    if not name.startswith("bioscrape."):
        return

    normalized = []
    previous_was_bullet = False
    for line in lines:
        stripped = line.strip()
        is_indented_detail = (
            (line[:1].isspace() and stripped)
            or stripped.startswith((":param", ":return", ":raises"))
        )
        if is_indented_detail:
            if normalized and normalized[-1].strip() and not previous_was_bullet:
                normalized.append("")
            normalized.append(f"* {stripped}")
            previous_was_bullet = True
            continue

        if previous_was_bullet and stripped:
            normalized.append("")
        normalized.append(line.rstrip())
        previous_was_bullet = False

    lines[:] = normalized


def setup(app):
    app.connect("autodoc-process-docstring", _normalize_legacy_docstrings)
