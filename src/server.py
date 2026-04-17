#!/usr/bin/env python3
"""MRBIGR2 aggregator — compatibility-mode single FastMCP exposing all 79 tools.

Imports each tool-mcp's ``register(mcp)`` and attaches them to one server.
Prefer running individual tool-mcps via ``mrbigr install <name>`` instead.
"""
import sys
import os
from pathlib import Path

_repo = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_repo / "src"))

for _mcp_dir in sorted((_repo / "tool-mcps").iterdir()):
    _srv = _mcp_dir / "src"
    if _srv.is_dir() and str(_srv) not in sys.path:
        sys.path.insert(0, str(_srv))

from fastmcp import FastMCP  # noqa: E402

HAS_FASTMCP = True
HAS_MODULES = True

_TOOL_MCPS = [
    "geno_mcp", "pheno_mcp", "gwas_mcp", "vis_mcp", "anno_mcp",
    "qtl_mcp", "net_mcp", "peak_mcp", "mr_mcp", "go_mcp",
]


def run_mcp_server() -> None:
    mcp = FastMCP("mrbigr2")

    for name in _TOOL_MCPS:
        srv_path = _repo / "tool-mcps" / name / "src" / "server.py"
        if not srv_path.exists():
            print(f"Warning: {name} server.py not found, skipping")
            continue
        import importlib.util
        spec = importlib.util.spec_from_file_location(f"{name}_server", str(srv_path))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        mod.register(mcp)

    mcp.run()


if __name__ == "__main__":
    run_mcp_server()
