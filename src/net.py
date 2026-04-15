#!/usr/bin/env python3
"""
NetMCP - Network Analysis Module
Pure Python implementation - No R dependencies

Based on MRBIGR/mrbigr/net.py

Functions:
- get_weight: Calculate edge weights from p-values
- module_identify: Network module detection (requires Java cluster_one)
- hub_identify: Hub gene identification
- module_network_plot: Network visualization

Dependencies:
- Java required for module_identify (cluster_one.jar)
- networkx for network analysis
"""
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from scipy.stats import expon
import subprocess
import os
import re
import warnings
import tempfile

# Try to import networkx for network analysis
try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    warnings.warn("networkx not installed. Network features limited.")

# Path to cluster_one jar
SCRIPT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLUSTER_ONE_JAR = os.path.join(SCRIPT_DIR, "utils", "cluster_one-1.0.jar")

# Check Java availability
def check_java():
    """Check if Java is available."""
    try:
        result = subprocess.run(['java', '-version'], capture_output=True, text=True, timeout=5)
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False

# Check cluster_one jar
def check_cluster_one():
    """Check if cluster_one jar exists."""
    return os.path.exists(CLUSTER_ONE_JAR)

HAS_JAVA = check_java()
HAS_CLUSTER_ONE = check_cluster_one()


def _normalize_edge_weight(edge_weight_input):
    if isinstance(edge_weight_input, pd.DataFrame):
        edge_weight = edge_weight_input.copy()
    else:
        edge_weight = pd.read_csv(edge_weight_input)

    rename_map = {}
    if 'source' in edge_weight.columns and 'row' not in edge_weight.columns:
        rename_map['source'] = 'row'
    if 'target' in edge_weight.columns and 'col' not in edge_weight.columns:
        rename_map['target'] = 'col'
    if rename_map:
        edge_weight = edge_weight.rename(columns=rename_map)

    required = {'row', 'col', 'weight'}
    if not required.issubset(edge_weight.columns):
        raise ValueError("edge_weight must contain row/col/weight columns")

    edge_weight = edge_weight[['row', 'col', 'weight']].copy()
    edge_weight['row'] = edge_weight['row'].astype(str)
    edge_weight['col'] = edge_weight['col'].astype(str)
    edge_weight['weight'] = pd.to_numeric(edge_weight['weight'], errors='coerce')
    edge_weight = edge_weight.dropna(subset=['row', 'col', 'weight'])
    edge_weight = edge_weight[edge_weight['weight'] > 0].reset_index(drop=True)
    return edge_weight


def _fallback_module_identify(edge_weight, module_size=3):
    if not HAS_NETWORKX:
        warnings.warn("networkx not installed. Module identification fallback unavailable.")
        return None

    if len(edge_weight) == 0:
        return pd.DataFrame(columns=['module', 'gene_num', 'genes'])

    G = nx.Graph()
    for _, row in edge_weight.iterrows():
        G.add_edge(row['row'], row['col'], weight=float(row['weight']))

    if G.number_of_edges() == 0:
        return pd.DataFrame(columns=['module', 'gene_num', 'genes'])

    try:
        communities = list(nx.algorithms.community.greedy_modularity_communities(G, weight='weight'))
    except Exception as e:
        warnings.warn(f"greedy_modularity_communities failed, falling back to connected components: {e}")
        communities = [set(c) for c in nx.connected_components(G)]

    modules = [sorted(list(c)) for c in communities if len(c) >= module_size]
    modules.sort(key=lambda genes: (-len(genes), genes[0] if genes else ""))

    return pd.DataFrame([
        {'module': idx + 1, 'gene_num': len(genes), 'genes': ' '.join(genes)}
        for idx, genes in enumerate(modules)
    ], columns=['module', 'gene_num', 'genes'])


# ========== Core Network Functions ==========

def get_weight(pvalue_list, pvalue_threshold=0.05):
    """Calculate edge weights from p-value list.
    
    Args:
        pvalue_list: DataFrame with mTrait, pTrait, pvalue columns
        pvalue_threshold: P-value threshold for significance
    
    Returns:
        DataFrame with row, col, weight (edge list format)
    """
    pvalue_list = pvalue_list[['mTrait', 'pTrait', 'pvalue']].copy()
    pvalue_list['pvalue'] = pd.to_numeric(pvalue_list['pvalue'], errors='coerce').clip(lower=1e-300, upper=1.0)
    pvalue_list = pvalue_list.dropna(subset=['mTrait', 'pTrait', 'pvalue'])
    
    # Create a symmetric p-value matrix across the union of all traits.
    id_list = sorted(set(pvalue_list['mTrait']) | set(pvalue_list['pTrait']))
    pvalue_matrix = pd.DataFrame(1.0, index=id_list, columns=id_list, dtype=float)

    pair_pvalues = (
        pvalue_list
        .groupby(['mTrait', 'pTrait'], as_index=False)['pvalue']
        .min()
    )
    for _, row in pair_pvalues.iterrows():
        m_trait = row['mTrait']
        p_trait = row['pTrait']
        pvalue = float(row['pvalue'])
        pvalue_matrix.loc[m_trait, p_trait] = min(pvalue_matrix.loc[m_trait, p_trait], pvalue)
        pvalue_matrix.loc[p_trait, m_trait] = min(pvalue_matrix.loc[p_trait, m_trait], pvalue)
    
    # Transform to weights
    pvalue_matrix = -np.log10(pvalue_matrix) + np.log10(pvalue_threshold)
    pvalue_matrix[pvalue_matrix < 0] = 0
    
    # Calculate weight using exponential distribution
    weight = 1 - expon.pdf(pvalue_matrix)
    weight = (weight + weight.T) / 2  # Make symmetric
    
    # Convert to edge list
    weight_csr = csr_matrix(np.triu(weight, k=1))
    weight_pair = []
    gene_id = pvalue_matrix.columns
    
    for row in range(weight.shape[0]):
        for col, w in zip(
            weight_csr.indices[weight_csr.indptr[row]:weight_csr.indptr[row + 1]],
            weight_csr.data[weight_csr.indptr[row]:weight_csr.indptr[row + 1]]
        ):
            weight_pair.append([gene_id[row], gene_id[col], w])
    
    return pd.DataFrame(weight_pair, columns=['row', 'col', 'weight'])


def module_identify(edge_weight_fn, module_size=3):
    """Identify network modules using cluster_one algorithm.
    
    Args:
        edge_weight_fn: Edge weight file path (CSV with row, col, weight)
        module_size: Minimum module size
    
    Returns:
        DataFrame with module, gene_num, genes columns
    """
    original_input_is_file = not isinstance(edge_weight_fn, pd.DataFrame)
    temp_file = None
    output_file = None

    try:
        edge_weight = _normalize_edge_weight(edge_weight_fn)

        if original_input_is_file:
            prefix = str(edge_weight_fn).replace('.csv', '').replace('.edge_list', '')
            output_file = f"{prefix}.cluster_one.result.csv"
        else:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
                edge_weight.to_csv(f, index=False)
                temp_file = f.name
            edge_weight_fn = temp_file

        if HAS_JAVA and HAS_CLUSTER_ONE:
            if output_file is None:
                prefix = edge_weight_fn.replace('.csv', '').replace('.edge_list', '')
                output_file = f"{prefix}.cluster_one.result.csv"

            cmd = f'java -jar {CLUSTER_ONE_JAR} -s {module_size} -f edge_list -F csv {edge_weight_fn} > {output_file} 2>/dev/null'
            subprocess.run(cmd, shell=True, check=True, capture_output=True, timeout=300)

            if os.path.exists(output_file):
                cluster_one_res = pd.read_csv(output_file)
                cluster_one_res = cluster_one_res.loc[cluster_one_res['P-value'] <= 0.05, :]
                cluster_one_res = cluster_one_res[['Size', 'Members']].sort_values(by='Size', ascending=False)
                cluster_one_res['module'] = np.arange(1, len(cluster_one_res) + 1)
                cluster_one_res.columns = ['gene_num', 'genes', 'module']
                return cluster_one_res[['module', 'gene_num', 'genes']]

        else:
            reason = "Java not available" if not HAS_JAVA else f"cluster_one.jar not found at {CLUSTER_ONE_JAR}"
            warnings.warn(f"{reason}. Falling back to NetworkX greedy modularity communities.")

        fallback_res = _fallback_module_identify(edge_weight, module_size=module_size)
        if fallback_res is not None and output_file and original_input_is_file:
            fallback_res.to_csv(output_file, index=False)
        return fallback_res

    except Exception as e:
        warnings.warn(f"Module identification failed: {e}")
        return None
    finally:
        if temp_file and os.path.exists(temp_file):
            os.unlink(temp_file)


def hub_identify(edge_weight, cluster_one_res):
    """Identify hub genes in network modules.
    
    Args:
        edge_weight: Edge weight DataFrame (row, col, weight)
        cluster_one_res: Module result from module_identify
    
    Returns:
        Tuple of (cluster_one_res with hub_gene, hub_res DataFrame)
    """
    if not HAS_NETWORKX:
        warnings.warn("networkx not installed. Hub identification skipped.")
        return cluster_one_res, pd.DataFrame()
    
    edge_weight = _normalize_edge_weight(edge_weight)
    cluster_one_res = cluster_one_res.copy()
    hub_res = pd.DataFrame()
    hub_list = []
    
    for _, row in cluster_one_res.iterrows():
        gene_list = row['genes'].split(' ')
        
        # Create subgraph
        m_ew = edge_weight[(edge_weight['row'].isin(gene_list)) & (edge_weight['col'].isin(gene_list))]
        
        if len(m_ew) == 0:
            hub_list.append('')
            continue
        
        # Build graph
        G = nx.Graph()
        for _, e in m_ew.iterrows():
            G.add_edge(e['row'], e['col'], weight=e['weight'])
        
        # Calculate hub scores
        try:
            hub_scores = nx.hits(G, max_iter=1000, normalized=True)[0]
        except Exception:
            degree_scores = nx.degree_centrality(G)
            max_score = max(degree_scores.values()) if degree_scores else 1.0
            hub_scores = {gene: (score / max_score if max_score else 0.0) for gene, score in degree_scores.items()}

        hub_df = pd.DataFrame({
            'gene_id': list(hub_scores.keys()),
            'hub_score': list(hub_scores.values())
        }).sort_values('hub_score', ascending=False)
        hub_res = pd.concat([hub_res, hub_df], ignore_index=True)

        top_hubs = hub_df[hub_df['hub_score'] >= 0.8]
        if top_hubs.empty and len(hub_df) > 0:
            top_hubs = hub_df.head(1)
        hub_str = ' '.join([f"{g}({s:.2f})" for g, s in zip(top_hubs['gene_id'], top_hubs['hub_score'])])
        hub_list.append(hub_str)
    
    cluster_one_res['hub_gene'] = hub_list
    return cluster_one_res, hub_res


def module_network_plot(edge_weight, cluster_one_res, hub_res=None, output_prefix="network", 
                       figsize=(12, 8), format='png'):
    """Plot network modules.
    
    Args:
        edge_weight: Edge weight DataFrame
        cluster_one_res: Module result DataFrame
        hub_res: Hub gene result DataFrame
        output_prefix: Output file prefix
        figsize: Figure size
        format: Output format (png, pdf)
    
    Returns:
        List of output file paths
    """
    if not HAS_NETWORKX:
        warnings.warn("networkx not installed. Network plotting skipped.")
        return []
    
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    output_files = []
    
    # Create hub score lookup
    hub_dict = {}
    if hub_res is not None and len(hub_res) > 0:
        for _, row in hub_res.iterrows():
            hub_dict[row['gene_id']] = row['hub_score']
    
    # Plot each module
    for idx, row in cluster_one_res.iterrows():
        gene_list = row['genes'].split(' ')
        
        # Get edges for this module
        m_ew = edge_weight[(edge_weight['row'].isin(gene_list)) & (edge_weight['col'].isin(gene_list))]
        
        if len(m_ew) == 0:
            continue
        
        # Create graph
        G = nx.Graph()
        for _, e in m_ew.iterrows():
            G.add_edge(e['row'], e['col'], weight=e['weight'])
        
        # Plot
        fig, ax = plt.subplots(figsize=figsize)
        
        # Node colors based on hub score
        node_colors = []
        node_sizes = []
        for node in G.nodes():
            score = hub_dict.get(node, 0.5)
            node_sizes.append(300 + score * 1000)
            if score >= 0.8:
                node_colors.append('red')
            elif score >= 0.3:
                node_colors.append('orange')
            else:
                node_colors.append('lightblue')
        
        # Draw network
        pos = nx.spring_layout(G, k=2, iterations=50)
        nx.draw_networkx_edges(G, pos, alpha=0.3, ax=ax)
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, ax=ax)
        nx.draw_networkx_labels(G, pos, font_size=8, ax=ax)
        
        ax.set_title(f"Module {row['module']}: {row['gene_num']} genes")
        ax.axis('off')
        
        # Save
        output_file = f"{output_prefix}_module{row['module']}.{format}"
        plt.savefig(output_file, bbox_inches='tight', dpi=150)
        plt.close()
        output_files.append(output_file)
    
    return output_files


def export_edge_list(weight_df, output_file):
    """Export edge list to file.
    
    Args:
        weight_df: DataFrame with row, col, weight
        output_file: Output file path
    
    Returns:
        Output file path
    """
    weight_df.to_csv(output_file, index=False)
    return output_file


def import_edge_list(edge_file):
    """Import edge list from file.
    
    Args:
        edge_file: Edge list CSV file
    
    Returns:
        DataFrame with row, col, weight
    """
    return pd.read_csv(edge_file)


def filter_edges_by_weight(edge_weight, min_weight=0.5):
    """Filter edges by minimum weight.
    
    Args:
        edge_weight: Edge weight DataFrame
        min_weight: Minimum weight threshold
    
    Returns:
        Filtered edge weight DataFrame
    """
    return edge_weight[edge_weight['weight'] >= min_weight].copy()


def get_network_stats(edge_weight):
    """Get network statistics.
    
    Args:
        edge_weight: Edge weight DataFrame with columns 'source', 'target', 'weight'
    
    Returns:
        Dictionary with network statistics
    """
    if not HAS_NETWORKX:
        return {'error': 'networkx not installed'}
    
    G = nx.Graph()
    for _, row in edge_weight.iterrows():
        source = row.get('source') or row.get('row') or row.get('from') or row.get(0)
        target = row.get('target') or row.get('col') or row.get('to') or row.get(1)
        weight = row.get('weight') or row.get(2) or 1.0
        G.add_edge(source, target, weight=weight)
    
    return {
        'n_nodes': G.number_of_nodes(),
        'n_edges': G.number_of_edges(),
        'density': nx.density(G),
        'is_connected': nx.is_connected(G),
        'n_components': nx.number_connected_components(G) if not nx.is_connected(G) else 1
    }


__all__ = [
    'get_weight', 'module_identify', 'hub_identify', 'module_network_plot',
    'export_edge_list', 'import_edge_list', 'filter_edges_by_weight', 'get_network_stats',
    'HAS_JAVA', 'HAS_CLUSTER_ONE', 'HAS_NETWORKX', 'check_java', 'check_cluster_one'
]
