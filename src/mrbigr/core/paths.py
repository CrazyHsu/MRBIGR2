"""Binary and root-path resolution for MRBIGR2.

Tool-mcps may be launched from arbitrary working directories, so every
call site that needs plink/gemma/FastTree/ClusterONE must go through
`find_binary()` rather than assuming cwd == repo root.

Resolution order:
  1. `MRBIGR_ROOT` environment variable, if set.
  2. Walk up from this file to find a directory containing `utils/` and
     `pyproject.toml`.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Optional


def repo_root() -> Path:
    """Return the MRBIGR2 repository root as a ``Path``."""
    env = os.environ.get("MRBIGR_ROOT")
    if env:
        return Path(env).resolve()

    here = Path(__file__).resolve()
    for candidate in (here, *here.parents):
        if (candidate / "utils").is_dir() and (candidate / "pyproject.toml").is_file():
            return candidate
    # Last resort: two levels up from src/mrbigr/core/paths.py == repo root
    return here.parents[3]


def utils_dir() -> Path:
    return repo_root() / "utils"


def find_binary(name: str) -> Optional[str]:
    """Locate a bundled or system binary by name.

    Looks first in ``<repo_root>/utils`` (where plink, gemma, FastTree,
    cluster_one-1.0.jar, and the Perl scripts live), then falls back to
    whatever ``$PATH`` provides.
    """
    bundled = utils_dir() / name
    if bundled.is_file() and os.access(bundled, os.X_OK):
        return str(bundled)
    return shutil.which(name)
