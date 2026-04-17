"""net_mcp — network analysis MCP (2 tools)."""
from __future__ import annotations

import sys
from pathlib import Path


def _bootstrap_syspath() -> None:
    here = Path(__file__).resolve()
    for candidate in here.parents:
        src_dir = candidate / "src" / "mrbigr"
        if src_dir.is_dir():
            src_parent = str(src_dir.parent)
            if src_parent not in sys.path:
                sys.path.insert(0, src_parent)
            return


_bootstrap_syspath()

from mrbigr.core import net  # noqa: E402


def register(mcp) -> None:  # type: ignore[no-untyped-def]
    @mcp.tool()
    def module_identify(edge_weight, module_size=3):
        """Identify network modules using ClusterONE or a NetworkX fallback."""
        import pandas as pd
        edge_df = pd.DataFrame(edge_weight) if isinstance(edge_weight, dict) else edge_weight
        result = net.module_identify(edge_df, module_size=module_size)
        return result.to_dict() if result is not None else None

    @mcp.tool()
    def hub_identify(edge_weight, cluster_one_res):
        """Identify hub genes in each network module."""
        import pandas as pd
        edge_df = pd.DataFrame(edge_weight) if isinstance(edge_weight, dict) else edge_weight
        cluster_df = pd.DataFrame(cluster_one_res) if isinstance(cluster_one_res, dict) else cluster_one_res
        cluster_res, hub_res = net.hub_identify(edge_df, cluster_df)
        return {
            "cluster_result": cluster_res.to_dict() if cluster_res is not None else None,
            "hub_result": hub_res.to_dict() if hub_res is not None else None,
        }


def main() -> None:
    from fastmcp import FastMCP
    mcp = FastMCP("net_mcp")
    register(mcp)
    mcp.run()


if __name__ == "__main__":
    main()
