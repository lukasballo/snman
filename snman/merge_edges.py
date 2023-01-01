import networkx as nx
from . import lanes
from . import graph_tools
from shapely import geometry, ops
from . import geometry_tools
import math
import numpy as np
import itertools as it


def merge_parallel_edges(street_graph):
    """
    Detect and merge all sets of edges sharing the same start/end nodes, incl. their attributes
    TODO: Merge this function with osmnx.utils_graph_get_digraph (including the merge of attributes)
    TODO: Avoid merging edges that are too far apart, e.g. parallel streets

    Parameters
    ----------
    street_graph : nx.MultiGraph

    Returns
    -------
    None
    """
    uv_index = {}
    for edge in street_graph.edges(data=True, keys=True):
        # Group edges by start/end nodes
        uv_key = ','.join([str(min(edge[0:2])), str(max(edge[0:2]))])
        if uv_key not in uv_index:
            uv_index[uv_key] = []
        uv_index[uv_key].append(edge)

    # Merge each set of grouped edges
    for uv_key, edge_list in uv_index.items():
        if len(edge_list) > 1:
            _merge_given_parallel_edges(street_graph, edge_list)
            pass


def _merge_given_parallel_edges(street_graph, edges):

    geometries = [edge[3].get('geometry') for edge in edges]
    offsets = np.array(geometry_tools._offset_distance(geometries))

    for index, edge in enumerate(edges):
        edge[3]['_geometry_offset'] = offsets[index]

    # sort the edges geometrically from left to the right
    edges_sorted_by_go = sorted(edges, key=lambda x: x[3].get('_geometry_offset'), reverse=True)

    # take the edge with the highest hierarchy as a parent edge
    i_parent_edge = 0
    edges_sorted_by_hierarchy = sorted(edges, key=lambda x: x[3].get('hierarchy'))
    parent_edge = edges_sorted_by_hierarchy[i_parent_edge]
    u = parent_edge[0]
    v = parent_edge[1]

    # merge the lanes from all
    parent_edge[3]['_merge_parallel_src_ln_desc'] = str([edge[3].get('ln_desc') for edge in edges_sorted_by_go])
    parent_edge[3]['ln_desc'] = list(it.chain(*[edge[3].get('ln_desc') for edge in edges_sorted_by_go]))

    # set the maximum maxspeed
    maxspeeds = [edge[3].get('maxspeed') for edge in edges]
    maxspeeds = list(filter(lambda e: isinstance(e, int), maxspeeds))
    maxspeed = max(maxspeeds) if len(maxspeeds) > 0 else None
    parent_edge[3]['maxspeed'] = maxspeed

    for index, edge in enumerate(edges_sorted_by_hierarchy):
        # remove edge, except if it's the parent edge
        if index != i_parent_edge:
            street_graph.remove_edges_from([edge])

    # Take the highest hierarchy level from all source edges
    parent_edge[3]['hierarchy'] = min([edge[3].get('hierarchy') for edge in edges])


def merge_consecutive_edges(street_graph):
    """
    Remove all intermediary nodes
    
    Parameters
    ----------
    street_graph : ox.MultiGraph
    
    Returns
    -------
    None
    """

    # Initialize 'consec_id' edge labels
    for edge in street_graph.edges(data=True, keys=True):
        edge[3]['consec_id'] = None

    # Label edges by consecutive clusters
    for node in street_graph.nodes:
        # Find all nodes having only two adjacent edges
        if street_graph.degree(node) == 2:
            edges = list(street_graph.edges(node, data=True, keys=True))
            # Start with the first adjacent edge and identify next edges to work on
            edge = edges[0]
            next_edges = _label_edge(street_graph, edge)
            # Repeat for each next_edge multiple times until entire clusters of consecutive edges are processed
            for i in range(3):
                next_edges = _label_edges(street_graph, next_edges)

    # Group edges by clusters
    edge_clusters = {}
    for edge in street_graph.edges(data=True, keys=True):
        consec_id = edge[3].get('consec_id')
        # Ignore edges that are in no cluster
        if consec_id is None:
            continue
        if consec_id not in edge_clusters:
            edge_clusters[consec_id] = []
        edge_clusters[consec_id].append(edge)

    # Merge edges within each cluster
    for key, edge_cluster in edge_clusters.items():
        _merge_given_consecutive_edges(street_graph, edge_cluster)
        pass


def _merge_given_consecutive_edges(street_graph, edges):

    nodes = []
    for edge in edges:
        nodes += edge[0:2]

    # Nodes that are adjacent to one edge -> outer nodes
    outer_nodes = [node for node in nodes if nodes.count(node) == 1]
    outer_nodes.sort()

    # Nodes that are adjacent to two edges -> middle nodes
    middle_nodes = [node for node in nodes if nodes.count(node) == 2]
    middle_nodes.sort()

    edge_chain = []

    # Stop here if the edge chain is corrupted
    # TODO: Why is this happening?
    if len(outer_nodes) != 2:
        #print('edge chain corrupted')
        return

    # Start with the first outer node
    node = outer_nodes[0]
    previous_edge = None

    # Do until reaching the other outer node
    while node != outer_nodes[1]:
        for edge in edges:

            # Ignore the previous edge
            if edge is previous_edge:
                continue

            # Append edge in ordinary direction
            elif edge[0] == node:
                edge_chain.append(edge)
                node = edge[1]

            # Append edge in reversed direction
            elif edge[1] == node:
                edge_chain.append(edge)
                node = edge[0]

    # Split the edge chain into subchains based on permitted modes
    edge_subchains = []
    motorized_access = None
    for edge in edge_chain:
        motorized_access_here = 'm>' in edge[3].get('ln_desc') or 'm<' in edge[3].get('ln_desc')\
                           or 'm-' in edge[3].get('ln_desc')
        if motorized_access != motorized_access_here:
            edge_subchains.append([])
        edge_subchains[-1].append(edge)
        motorized_access = motorized_access_here

    # Merge all edges within each subchain
    for edge_subchain in edge_subchains:

        if len(edge_subchain) < 2:
            continue

        # Find subchain outer nodes
        subchain_nodes = []
        for edge in edge_subchain:
            subchain_nodes += edge[0:2]
        subchain_outer_nodes = list(set([node for node in subchain_nodes if subchain_nodes.count(node) == 1]))
        subchain_outer_nodes.sort()
        subchain_middle_nodes = list(set([node for node in subchain_nodes if subchain_nodes.count(node) > 1]))
        subchain_middle_nodes.sort()

        # Reverse edge geometries if necessary
        node = subchain_outer_nodes[0]
        for edge in edge_subchain:
            edge[3]['__previous_node'] = node
            if node == edge[0]:
                node = edge[1]
            elif node == edge[1]:
                graph_tools._reverse_edge(street_graph, edge, reverse_topology=False)
                node = edge[0]
            edge[3]['__next_node'] = node

        # Find the longest edge
        sorted_edges = sorted(edge_subchain, key=lambda x: x[3].get('length'))
        longest_edge = sorted_edges[-1]

        # Initialize the merged edge based on the longest edge
        merged_data = longest_edge[3].copy()

        # Merge geometries
        geometries = [edge[3]['geometry'] for edge in edge_subchain]
        for i in range(2):
            geometries = [
                (lambda geom: (ops.linemerge(geom) if geom.geom_type=='MultiLineString' else geom))(geom)
                for geom in geometries
            ]
        multi_line = geometry.MultiLineString(geometries)
        merged_line = ops.linemerge(multi_line)
        merged_data['geometry'] = merged_line
        # Update length
        merged_data['length'] = merged_line.length
        merged_data['__outer_nodes'] = str(subchain_outer_nodes)
        merged_data['__middle_nodes'] = str(subchain_middle_nodes)
        merged_data['__nodes'] = str(subchain_nodes)
        merged_data['__n_edges'] = str(len(edge_subchain))

        # Delete the old edges
        for edge in edge_subchain:
            street_graph.remove_edge(edge[0], edge[1], edge[2])
            pass

        # Delete the intermediary nodes
        for node in subchain_middle_nodes:
            street_graph.remove_node(node)

        # Create a new merged edge
        street_graph.add_edge(subchain_outer_nodes[0], subchain_outer_nodes[1], **merged_data)


def _label_edge(graph, edge):
    edge_data = edge[3]
    neighbors_groups = graph_tools._get_neighbors(graph, edge)
    edge[3]['neighbors'] = str([len(group) for group in neighbors_groups])
    neighbors_groups.append([edge])
    next_edges = []
    consecutive_cluster_id = None

    for neighbors_group in neighbors_groups:
        # Discard neighbor groups with multiple edges, indicating that the node has >2 adjacent edges
        if len(neighbors_group) != 1:
            continue

        neighbor = neighbors_group[0]
        # Add to next_edges if it has not been processed yet (i.e. has no consec_id)
        neighbor[3].get('consec_id') or next_edges.append(neighbor)
        # If it has a consec_id, take it over
        consecutive_cluster_id = consecutive_cluster_id or neighbor[3].get('consec_id')

    # If no consec_id was found, assign a new one based on the own edge_hash
    consecutive_cluster_id = consecutive_cluster_id or graph_tools._edge_hash(edge)

    # Write the consec_id into all edges of the cluster
    for edge in next_edges:
        #print(consecutive_cluster_id)
        edge[3]['consec_id'] = consecutive_cluster_id
        pass

    return next_edges


def _label_edges(graph, edges):
    next_edges = []
    for edge in edges:
        result = _label_edge(graph, edge)
        if (len(result) > 0):
            next_edges += result

    return next_edges


def resolve_one_sided_intersections(graph):
    """
    Connect one-sided intersections on roads with multiple separate edges to all edges
    """

    for node in graph.nodes(data=True, keys=True):
        pass
