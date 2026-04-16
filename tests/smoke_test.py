"""Phase-gating smoke tests for the MRBIGR2 refactor.

Each phase of feat/real-mcp-refactor must keep these tests green. They exercise
the boot path of the FastMCP server plus a read-only import of every domain
module, without touching data/ files.

Run: pytest tests/smoke_test.py -q
"""
from __future__ import annotations

import importlib
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC = REPO_ROOT / "src"

# Domain modules that must stay importable through every phase.
# In Phase 1 they live at src/<name>.py; in Phase 2+ they move to
# src/mrbigr/core/<name>.py and this fixture will be updated in lockstep.
DOMAIN_MODULES = [
    "pheno", "geno", "gwas", "vis", "anno",
    "qtl", "mr", "go", "net", "peak", "multi",
]


def _import_from(path: Path, name: str):
    sys.path.insert(0, str(path))
    try:
        return importlib.import_module(name)
    finally:
        sys.path.remove(str(path))


@pytest.mark.parametrize("module", DOMAIN_MODULES)
def test_domain_module_imports(module: str) -> None:
    """Every domain module must import without side effects."""
    mod = _import_from(SRC, module)
    assert mod is not None


def test_server_boots_and_exits() -> None:
    """Run `python src/server.py --help`-equivalent: import the module and
    confirm FastMCP and all domain imports succeed. We do *not* actually
    enter the stdio loop.
    """
    result = subprocess.run(
        [sys.executable, "-c",
         "import sys; sys.path.insert(0, 'src'); "
         "import server; "
         "assert server.HAS_FASTMCP, 'fastmcp missing'; "
         "assert server.HAS_MODULES, 'domain modules missing'; "
         "print('OK')"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, (
        f"server.py failed to import.\nSTDOUT: {result.stdout}\n"
        f"STDERR: {result.stderr}"
    )
    assert "OK" in result.stdout
