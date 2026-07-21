from __future__ import annotations

import importlib
from importlib import metadata as importlib_metadata
import inspect
import re
import sys
import sysconfig
from pathlib import Path

import sphinx

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

import bioscrape

project = "Bioscrape"
author = "Biocircuits, California Institute of Technology"
copyright = "2026, Biocircuits, California Institute of Technology"

try:
    release = importlib_metadata.version("Bioscrape")
except importlib_metadata.PackageNotFoundError:
    release = "1.3.0"

version = ".".join(release.split(".", 2)[:2])

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.linkcode",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx_copybutton",
    "sphinx_math_dollar",
    "sphinx_toggleprompt",
    "nbsphinx",
    "nbsphinx_link",
    "recommonmark",
    "numpydoc",
]

source_suffix = [".rst"]

autosummary_generate = True
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
    "exclude-members": "__init__, __weakref__, __repr__, __str__, __hash__",
}

# For classes, include the class docstring instead of pulling in init details.
autoclass_content = "class"

autodoc_mock_imports = [
    "roadrunner",
]

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "**.ipynb_checkpoints",
    "Thumbs.db",
    ".DS_Store",
    "user_guide/model_api.rst",
    "user_guide/bioscrape_xml.rst",
]

pygments_style = "sphinx"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}
intersphinx_disabled_reftypes = ["py:keyword"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_title = "Bioscrape documentation"

# Match BioCRNpyler's convention: single backticks are Python objects and will
# link into the API reference when Sphinx can resolve them.
default_role = "py:obj"

sphinx_version = tuple(int(x) for x in sphinx.__version__.split(".")[:2])
if sphinx_version >= (4, 0):
    mathjax3_config = {
        "tex": {
            "inlineMath": [["\\(", "\\)"]],
            "displayMath": [["\\[", "\\]"]],
        }
    }
else:
    mathjax_config = {
        "tex2jax": {
            "inlineMath": [["\\(", "\\)"]],
            "displayMath": [["\\[", "\\]"]],
        },
    }

copybutton_prompt_text = r">>> |\.\.\. "
copybutton_prompt_is_regexp = True

nbsphinx_execute = "never"
nbsphinx_allow_errors = True
nbsphinx_kernel_name = "python3"

# The examples are tutorial artifacts; their output cells are preserved and
# source links are more useful than long execution traces in RTD builds.
nbsphinx_timeout = 120

numpydoc_show_class_members = False
numpydoc_class_members_toctree = False

napoleon_use_ivar = False
napoleon_use_param = True
napoleon_preprocess_types = True
napoleon_custom_sections = [
    ("Attributes", "params_style"),
]

_TYPE_ALIAS_MODULES = [
    "bioscrape.types",
    "bioscrape.simulator",
    "bioscrape.analysis",
    "bioscrape.random",
    "bioscrape.inference",
    "bioscrape.inference_setup",
    "bioscrape.pid_interfaces",
    "bioscrape.sbmlutil",
]


def _build_type_aliases() -> dict[str, str]:
    aliases: dict[str, str] = {}
    for module_name in _TYPE_ALIAS_MODULES:
        try:
            module = importlib.import_module(module_name)
        except Exception:
            continue

        for class_name, obj in inspect.getmembers(module, inspect.isclass):
            obj_module = getattr(obj, "__module__", module_name)
            if obj_module.startswith("bioscrape"):
                aliases.setdefault(
                    class_name, f":class:`~{obj_module}.{class_name}`"
                )
    return aliases


napoleon_type_aliases = _build_type_aliases()
autodoc_type_aliases = {
    name: napoleon_type_aliases[name] for name in sorted(napoleon_type_aliases)
}

#
# Docstring pre-processing
#

eqn_substitutions = [
    (r"<-->", r"\\rightleftharpoons"),
    (r"-->", r"\\rightarrow"),
    (r"\.\.\.", r"\\dots"),
    (r":", r"\\mathord{:}"),
    (r"\{\}", r"\\emptyset"),
    (r" >> ", r" \\gg "),
    (r" << ", r" \\ll "),
    (r"'([\w -]+)'", r"{\\text{\1}}"),
    (r"\[([\w -]+)\]", r"[\\text{\1}]"),
    (r"^[ ]+", r""),
    (r"\$[ ]+", r"$"),
    (r"&    ", r"& \\qquad"),
]

txt_substitutions = [
    (r"<-->", r"$\\rightleftharpoons$"),
    (r"-->", r"$\\rightarrow$"),
    (r"\{\}", r"$\\emptyset$"),
    (r" >> ", r" $\\gg$ "),
    (r" << ", r" $\\ll$ "),
]


def _process_string(s: str, substitutions: list[tuple[str, str]]) -> str:
    for pattern, repl in substitutions:
        s = re.sub(pattern, repl, s)
    return s


def _has_numpydoc_sections(lines) -> bool:
    section_names = {
        "Attributes",
        "Examples",
        "Notes",
        "Parameters",
        "Raises",
        "Returns",
        "See Also",
        "Warnings",
        "Yields",
    }
    for index, line in enumerate(lines[:-1]):
        underline = lines[index + 1].strip()
        if line.strip() in section_names and underline and set(underline) == {"-"}:
            return True
    return False


def _normalize_legacy_docstrings(app, what, name, obj, options, lines):
    """Make older Cython docstrings acceptable to docutils."""
    if not name.startswith("bioscrape.") or _has_numpydoc_sections(lines):
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


def _preprocess_docstrings(app, what, name, obj, options, lines):
    """Render lightweight reaction notation in docstrings as math."""
    in_equation = False
    for i, line in enumerate(lines):
        eqn_list = list(re.finditer(r"\$[^$]+\$", line))
        if re.match(r"^[ ]*\$\$$", line):
            in_equation = not in_equation

        if in_equation:
            line = _process_string(line, eqn_substitutions)
        elif eqn_list:
            line, offset = "", 0
            for match in eqn_list:
                line += _process_string(
                    lines[i][offset : match.start()], txt_substitutions
                )
                offset = match.end()
                line += _process_string(match.group(0), eqn_substitutions)

            line += _process_string(
                lines[i][eqn_list[-1].end() :], txt_substitutions
            )
        else:
            line = _process_string(line, txt_substitutions)

        lines[i] = line


# -----------------------------------------------------------------------------
# Source code links
# -----------------------------------------------------------------------------


def _source_path_for_binary(path: Path, modname: str) -> Path | None:
    module_name = modname.rsplit(".", 1)[-1]
    binary_suffixes = {".pyd", ".so", ".dll", ".dylib"}
    if extension_suffix:
        binary_suffixes.add(extension_suffix)

    if not any(path.name.endswith(suffix) for suffix in binary_suffixes):
        return None

    for source_suffix in (".pyx", ".py"):
        candidate = ROOT / "bioscrape" / f"{module_name}{source_suffix}"
        if candidate.exists():
            return candidate
    return None


def linkcode_resolve(domain, info):
    """Resolve documented Python objects to their GitHub source location."""
    if domain != "py":
        return None

    modname = info.get("module")
    fullname = info.get("fullname")
    if not modname or not fullname:
        return None

    submod = sys.modules.get(modname)
    if submod is None:
        try:
            submod = importlib.import_module(modname)
        except Exception:
            return None

    obj = submod
    for part in fullname.split("."):
        try:
            obj = getattr(obj, part)
        except Exception:
            return None

    try:
        obj = inspect.unwrap(obj)
    except Exception:
        pass

    try:
        filename = inspect.getsourcefile(obj) or inspect.getfile(obj)
    except Exception:
        return None

    if not filename:
        return None

    module = inspect.getmodule(obj)
    if module is not None and not module.__name__.startswith("bioscrape"):
        return None

    try:
        source, lineno = inspect.getsourcelines(obj)
    except Exception:
        source, lineno = None, None

    try:
        path = Path(filename).resolve()
    except OSError:
        return None

    source_path = _source_path_for_binary(path, modname)
    if source_path is not None:
        path = source_path.resolve()
        source, lineno = None, None

    relpath = None
    for base in (ROOT, Path(bioscrape.__file__).resolve().parents[1]):
        try:
            relpath = path.relative_to(base)
            break
        except ValueError:
            continue

    if relpath is None:
        parts = path.parts
        if "bioscrape" not in parts:
            return None
        relpath = Path(*parts[parts.index("bioscrape") :])

    if not relpath.parts or relpath.parts[0] != "bioscrape":
        return None

    linespec = ""
    if lineno and source:
        linespec = f"#L{lineno}-L{lineno + len(source) - 1}"

    return (
        "https://github.com/biocircuits/bioscrape/blob/main/"
        f"{relpath.as_posix()}{linespec}"
    )


doctest_global_setup = """
import numpy as np
"""


def setup(app):
    app.connect(
        "autodoc-process-docstring", _normalize_legacy_docstrings, priority=400
    )
    app.connect("autodoc-process-docstring", _preprocess_docstrings, priority=410)
