"""Load mcps.yaml and drive the install/register lifecycle.

A 5-minute in-process cache prevents redundant disk + `claude mcp list`
round-trips during multi-MCP operations (e.g. installing everything in
one shell session).
"""
from __future__ import annotations

import time
from pathlib import Path
from typing import Any

import yaml

from ..core.paths import repo_root
from .mcp import MCP, Status

CACHE_TTL_SECONDS = 300

_ALLOWED_FIELDS = set(MCP.__dataclass_fields__.keys()) - {"name"}


class MCPManager:
    """Read-through manager for the MCP registry."""

    def __init__(self, root: Path | None = None) -> None:
        self.root = (root or repo_root()).resolve()
        self.registry_path = self.root / "mcps.yaml"
        self._cache: dict[str, MCP] | None = None
        self._cache_ts: float = 0.0

    # ---- loading -------------------------------------------------------

    def _fresh_load(self) -> dict[str, MCP]:
        if not self.registry_path.is_file():
            return {}
        with open(self.registry_path) as f:
            data: dict[str, Any] = yaml.safe_load(f) or {}
        mcps: dict[str, MCP] = {}
        for name, entry in (data.get("mcps") or {}).items():
            fields = {k: v for k, v in (entry or {}).items() if k in _ALLOWED_FIELDS}
            mcps[name] = MCP(name=name, **fields)
        return mcps

    def load(self, force: bool = False) -> dict[str, MCP]:
        now = time.time()
        if not force and self._cache is not None and (now - self._cache_ts) < CACHE_TTL_SECONDS:
            return self._cache
        self._cache = self._fresh_load()
        self._cache_ts = now
        return self._cache

    def invalidate(self) -> None:
        self._cache = None

    # ---- queries -------------------------------------------------------

    def get(self, name: str) -> MCP:
        mcps = self.load()
        if name not in mcps:
            raise KeyError(f"Unknown MCP: {name}. Known: {sorted(mcps)}")
        return mcps[name]

    def list(self) -> list[MCP]:
        return list(self.load().values())

    def status(self, name: str) -> Status:
        return self.get(name).status(self.root)

    # ---- actions -------------------------------------------------------

    def install(self, name: str, python_bin: str | None = None) -> Status:
        mcp = self.get(name)
        mcp.register(self.root, python_bin=python_bin)
        self.invalidate()
        return mcp.status(self.root)

    def uninstall(self, name: str) -> Status:
        mcp = self.get(name)
        mcp.unregister()
        self.invalidate()
        return mcp.status(self.root)
