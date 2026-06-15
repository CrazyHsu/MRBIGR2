"""Helpers for MCP arg robustness.

The MCP boundary can deliver structured (``dict``/``list``) arguments to
untyped tool parameters as JSON *strings*. ``maybe_json_loads`` decodes such a
string back to a Python object so downstream coercion (``pd.DataFrame`` /
``np.array`` / ``.columns`` ...) works. Real file-path strings never start with
``[`` or ``{`` and are therefore returned unchanged so path-based loaders still
apply.
"""
import gzip
import json


def maybe_json_loads(value):
    """If ``value`` is a JSON object/array *string*, decode it; else return as-is."""
    if isinstance(value, str):
        s = value.strip()
        if s[:1] in "[{":
            try:
                return json.loads(s)
            except Exception:
                return value
    return value


def coerce_table(value):
    """Coerce an MCP arg (json-string / dict / list / CSV path / DataFrame) to a DataFrame."""
    import os
    import pandas as pd

    value = maybe_json_loads(value)
    if isinstance(value, pd.DataFrame):
        return value
    if isinstance(value, str):
        if os.path.isfile(value):
            return pd.read_csv(value, sep=None, engine="python")
        return value
    return pd.DataFrame(value)


def open_text_maybe_gz(path):
    """Open a text file, transparently handling gzip (.gz / gzip magic)."""
    p = str(path)
    is_gz = p.endswith(".gz")
    if not is_gz:
        try:
            with open(p, "rb") as fh:
                is_gz = fh.read(2) == b"\x1f\x8b"
        except OSError:
            is_gz = False
    return gzip.open(p, "rt") if is_gz else open(p)
