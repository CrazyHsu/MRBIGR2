"""Shared utilities for the mrbigr meta-orchestrator."""
from __future__ import annotations

import os
import sys


def ensure_repo_on_syspath() -> None:
    """Legacy-compat: for direct `python src/server.py` runs without `pip install -e .`.

    After `pip install -e .`, this is a no-op — but during development
    people do still invoke the server as a bare script, and the mrbigr.*
    package would otherwise be unreachable.
    """
    from .core.paths import repo_root
    src = str(repo_root() / "src")
    if src not in sys.path:
        sys.path.insert(0, src)
