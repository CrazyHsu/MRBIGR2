"""Input coercion tests for the GO MCP wrapper.

Regression for plot_gsea_results: an inline GSEA result delivered as a JSON
*string* (how this MCP client passes untyped object/array args) must be decoded
via coerce_table, not crash with ``'str' object has no attribute 'copy'``.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SERVER_PY = REPO_ROOT / "tool-mcps" / "go_mcp" / "src" / "server.py"
REPO_SRC = REPO_ROOT / "src"
if str(REPO_SRC) not in sys.path:
    sys.path.insert(0, str(REPO_SRC))


def _load_go_server():
    spec = importlib.util.spec_from_file_location("go_mcp_test_server", SERVER_PY)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class FakeMCP:
    def __init__(self) -> None:
        self.tools = {}

    def tool(self):
        def decorate(fn):
            self.tools[fn.__name__] = fn
            return fn

        return decorate


def test_plot_gsea_results_accepts_json_string(tmp_path: Path) -> None:
    server = _load_go_server()
    fake_mcp = FakeMCP()
    server.register(fake_mcp)

    # The inline list arrives as a JSON string at the untyped tool parameter.
    gsea_json = (
        '[{"Term":"t1","NES":2.0,"FDR q-val":0.01,"NOM p-val":0.001},'
        '{"Term":"t2","NES":-1.5,"FDR q-val":0.04,"NOM p-val":0.01}]'
    )

    out = fake_mcp.tools["plot_gsea_results"](
        gsea_json, str(tmp_path / "gsea"), plot_mode="nes_barplot"
    )

    assert isinstance(out, dict)
    assert out  # at least one plot produced
    assert any(Path(p).is_file() for p in out.values())
