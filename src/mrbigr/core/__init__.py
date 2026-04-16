"""Shared analytical core for MRBIGR2.

Every tool-mcp under tool-mcps/ and the legacy aggregator at src/server.py
imports from this package. Modules here are domain code; they must not
depend on FastMCP or any agent-orchestration concern.
"""

from . import (
    anno,
    geno,
    go,
    gwas,
    mr,
    net,
    parallel,
    peak,
    pheno,
    qtl,
    vis,
)

__all__ = [
    "anno", "geno", "go", "gwas", "mr", "net",
    "parallel", "peak", "pheno", "qtl", "vis",
]
