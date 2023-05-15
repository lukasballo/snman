import networkx as nx
from . import graph_utils, geometry_tools, utils, street_graph, space_allocation
from. constants import *
from shapely import geometry, ops
import numpy as np
import itertools as it


def merge_parallel_edges(G):
    """
    Detect and merge all sets of edges sharing the same start/end nodes, incl. their attributes
    TODO: Avoid merging edges that are too far apart, e.g. parallel streets

    Parameters
    ----------
    G : nx.MultiDiGraph
        street graph

    Returns
    -------
    None
    """

    edge_lists = {}
    for edge in G.edges(data=True, keys=True):
        # group edges by start/end nodes and layer (uvl = u, v, layer)
        uvl = (min(edge[0:2]), max(edge[0:2]), edge[3].get('layer'))
        # initialize empty edge list if necessary
        if uvl not in edge_lists:
            edge_lists[uvl] = []
        edge_lists[uvl].append(edge)

    # Merge each set of grouped edges
    for uvl, edge_list in edge_lists.items():
        if len(edge_list) > 1:
            _merge_given_parallel_edges(G, *uvl, edge_list)
            pass


def _merge_given_parallel_edges(G, u, v, l, edges):
    """
    Merge the parallel edges given in a list

    Parameters
    ----------
    G : nx.MultiDiGraph
        street graph
    edges : list
        list of edges to be merged

    Returns
    -------
    None
    """

    # normalize edge directions
    edges_normalized = []
    for edge in edges:
        if edge[0:2] == (v, u):
            reversed_edge = street_graph.reverse_edge(G, *edge[0:3])
            edges_normalized.append(reversed_edge)
        else:
            edges_normalized.append(edge)

    edges = edges_normalized

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

    # merge the lanes from all
    parent_edge_data = G.edges[parent_edge[0:3]]
    parent_edge_data['_merge_parallel_src_ln_desc'] = str([edge[3].get('ln_desc') for edge in edges_sorted_by_go])
    parent_edge_data['ln_desc'] = list(it.chain(*[edge[3].get('ln_desc') for edge in edges_sorted_by_go]))

    # merge sensors
    parent_edge_data['sensors_forward'] = list(set(
        utils.flatten_list(
            [edge[3].get('sensors_forward', []) for edge in edges_sorted_by_go]
        )
    ))
    parent_edge_data['sensors_backward'] = list(set(
        utils.flatten_list(
            [edge[3].get('sensors_backward', []) for edge in edges_sorted_by_go]
        )
    ))

    # set the maximum maxspeed
    maxspeeds = [edge[3].get('maxspeed') for edge in edges]
    maxspeeds = list(filter(lambda e: isinstance(e, int), maxspeeds))
    maxspeed = max(maxspeeds) if len(maxspeeds) > 0 else None
    parent_edge_data['maxspeed'] = maxspeed

    for index, edge in enumerate(edges_sorted_by_hierarchy):
        # remove edge, except if it's the parent edge
        if index != i_parent_edge:
            G.remove_edges_from([edge])

    # Take the highest hierarchy level from all source edges
    parent_edge_data['hierarchy'] = min([edge[3].get('hierarchy') for edge in edges])


def merge_consecutive_edges(G):

    # create a subgraph that only contains degree=2 nodes
    H = G.copy()
    degrees = dict(H.degree)
    for node in list(H.nodes()):
        if degrees[node] != 2:
            H.remove_node(node)

    # identify weakly connected components as clusters to be merged
    wccs = list(nx.weakly_connected_components(H))
    clusters = {}
    for i, wcc in enumerate(wccs):
        for node in wcc:
            H.nodes[node]['_cc_id'] = i
            # initialize a list for edges
            clusters[i] = set()

    # write the cluster ids to the main graph
    for i, data in H.nodes.items():
        cc_id = data['_cc_id']
        edges = set(G.in_edges(i, keys=True)).union(set(G.out_edges(i, keys=True)))
        clusters[cc_id].update(edges)
        for edge in edges:
            nx.set_edge_attributes(G, {edge: {'_cc_id': cc_id}})

    # bring the edges in each cluster into a consistent direction
    for i, edges in clusters.items():
        # make a small graph for each cluster, it must be undirected so that there is a shortest path
        H = nx.MultiGraph(edges)
        outer_nodes = [node for node, degree in H.degree() if degree == 1]
        if len(outer_nodes) != 2:
            continue
        sp = nx.shortest_path(H, source=outer_nodes[0], target=outer_nodes[1])
        # reconstruct the path edges from the shortest path nodes
        path_pairs = list(zip(sp[:-1], sp[1:]))
        path_keys = [list(H[pair[0]][pair[1]].keys())[0] for pair in path_pairs]
        path_edges = list(zip(sp[:-1], sp[1:], path_keys))
        # reverse those edges that point in one direction and leave the others as they are
        edges_to_merge = []
        for edge in path_edges:
            if edge not in G.edges:
                edges_to_merge.append(
                    street_graph.reverse_edge(G, edge[1], edge[0], edge[2])
                )
            else:
                edges_to_merge.append(
                    (*edge, G.edges[edge[0:3]])
                )

        # merge the edges
        _merge_given_consecutive_edges(G, edges_to_merge)


def _merge_given_consecutive_edges(G, edge_chain):
    """
    Merge the consecutive edges in a given list

    Parameters
    ----------
    G : nx.MultiGraph
        street graph
    edges : list
        list of edges to be merged

    Returns
    -------
    None
    """

    nodes = []
    for edge in edge_chain:
        nodes += edge[0:2]

    # split the edge chain into subchains based on permitted modes
    edge_subchains = []
    motorized_access = None
    for edge in edge_chain:
        ls = space_allocation._lane_stats(KEY_LANES_DESCRIPTION)
        motorized_access_here = 0 < (
            ls.n_lanes_motorized_forward
            + ls.n_lanes_motorized_backward
            + ls.n_lanes_motorized_both_ways
            + ls.n_lanes_motorized_direction_tbd
        )
        if motorized_access != motorized_access_here:
            edge_subchains.append([])
        edge[3]['_subchain'] = len(edge_subchains)
        edge_subchains[-1].append(edge)
        motorized_access = motorized_access_here

    # merge all edges within each subchain
    for edge_subchain in edge_subchains:

        if len(edge_subchain) < 2:
            continue

        # Find the longest edge
        sorted_edges = sorted(edge_subchain, key=lambda x: x[3].get('length', 0))
        longest_edge = sorted_edges[-1]

        # Initialize the merged edge based on the longest edge
        merged_data = longest_edge[3].copy()

        # Merge geometries
        geometries = [edge[3]['geometry'] for edge in edge_subchain]
        for i in range(2):
            geometries = [
                (lambda geom: (ops.linemerge(geom) if geom.geom_type == 'MultiLineString' else geom))(geom)
                for geom in geometries
            ]
        multi_line = geometry.MultiLineString(geometries)
        merged_line = ops.linemerge(multi_line)
        merged_data['geometry'] = merged_line
        # Update length
        merged_data['length'] = merged_line.length

        # Merge other attributes
        merged_data['sensors_forward'] = list(set(
            utils.flatten_list(
                [edge[3].get('sensors_forward', []) for edge in edge_subchain]
            )
        ))
        merged_data['sensors_backward'] = list(set(
            utils.flatten_list(
                [edge[3].get('sensors_backward', []) for edge in edge_subchain]
            )
        ))

        # Delete the old edges
        for edge in edge_subchain:
            G.remove_edge(*edge[0:3])

        # Create a new merged edge
        u = edge_subchain[0][0]
        v = edge_subchain[-1][1]
        G.add_edge(u, v, **merged_data)
