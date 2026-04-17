# net_mcp

Network analysis MCP — 2 tools backed by `mrbigr.core.net`.

## Tools

`module_identify` (ClusterONE with NetworkX fallback), `hub_identify`.

Needs Java + `utils/cluster_one-1.0.jar` for the canonical ClusterONE
path; NetworkX-only fallback is used when Java is absent.

Standalone: `python tool-mcps/net_mcp/src/server.py`
Register:   `mrbigr install net_mcp`
