import networkx as nx
from . import graph, geometry_tools, utils, street_graph, space_allocation
from. constants import *
from shapely import geometry, ops
import shapely
import numpy as np
import itertools as it
import copy


def merge_parallel_edges(G, max_hausdorff_distance=50):
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
        data = edge[3]
        # skip this edge if it should not be included in simplification
        if not data.get('_include_in_simplification', True):
            continue
        # group edges by start/end nodes and layer (uvl = u, v, layer)
        uvl = (min(edge[0:2]), max(edge[0:2]), edge[3].get('layer'))
        # initialize empty edge list if necessary
        if uvl not in edge_lists:
            edge_lists[uvl] = []
        edge_lists[uvl].append(edge)

    # Merge each set of grouped edges
    for uvl, edge_list in edge_lists.items():
        if len(edge_list) > 1:
            _merge_given_parallel_edges(G, *uvl, edge_list, max_hausdorff_distance)
            pass


def _merge_given_parallel_edges(G, u, v, l, edges, max_hausdorff_distance):
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

    # take the edge with the highest osm highway tag level
    i_parent_edge = 0
    edges_sorted_by_hierarchy = sorted(edges, key=lambda x: OSM_HIGHWAY_VALUES[x[3]['highway']]['level'])
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
            graph.safe_remove_edge(G, *edge[0:3])


def merge_consecutive_edges(G, distinction_attributes=set()):

    # create a subgraph that only contains nodes that are all of these
    # - degree=2
    # - included in the simplification
    # - have only one layer
    H = G.copy()
    degrees = dict(H.degree)
    for n, data in G.nodes.items():
        if degrees[n] == 2 and data.get('_include_in_simplification', True) and len(data.get('layers', [])) == 1:
            pass
        else:
            H.remove_node(n)

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
        _merge_given_consecutive_edges(G, edges_to_merge, distinction_attributes)


def _merge_given_consecutive_edges(G, edge_chain, distinction_attributes):
    """
    Merge the consecutive edges in a given list

    Parameters
    ----------
    G : nx.MultiGraph
        street graph
    edge_chain : list
        edges to be merged
    distinction_attributes : set
        if these attributes differ, then start a new subchain,
        e.g., when the lanes are different

    Returns
    -------
    None
    """

    nodes = []
    for edge in edge_chain:
        nodes += edge[0:2]

    # split the edge chain into subchains based on permitted modes
    edge_subchains = []
    # initialize distinction variables
    motorized_access = None
    distinction_attr_values = {attr: None for attr in distinction_attributes}
    for edge in edge_chain:
        ls = space_allocation._lane_stats(edge[3].get(KEY_LANES_DESCRIPTION, []))
        motorized_access_here = 0 < (
            ls.n_lanes_motorized_forward
            + ls.n_lanes_motorized_backward
            + ls.n_lanes_motorized_both_ways
            + ls.n_lanes_motorized_direction_tbd
        )
        distinction_attr_values_here = {attr: edge[3][attr] for attr in distinction_attributes}
        # start a new subchain if the distinction triggers
        if motorized_access != motorized_access_here or distinction_attr_values != distinction_attr_values_here:
            edge_subchains.append([])
        edge[3]['_subchain'] = len(edge_subchains)
        # add the edge to the last subchain
        edge_subchains[-1].append(edge)
        motorized_access = motorized_access_here
        distinction_attr_values = distinction_attr_values_here

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

        # add a list of merged intermediary nodes
        merged_data['_intermediary_nodes'] = [edge[0] for edge in edge_subchain[1:]]

        # Create a new merged edge
        u = edge_subchain[0][0]
        v = edge_subchain[-1][1]
        G.add_edge(u, v, **merged_data)

        # Delete the old edges
        for edge in edge_subchain:
            graph.safe_remove_edge(G, *edge[0:3])


def reconstruct_consecutive_edges(G):
    """
    Converts all merged edges into their original parts.
    This function is needed to write back the rebuilt lanes into the original street graph if we merge
    consecutive edges in the subgraph (e.g., when rebuilding only the main roads).

    For simplicity, the reconstructed graph does not maintain the original locations of the intermediary points,
    nor the real edge geometries. Instead, the intermediary points are distributed evenly and connected with straight
    lines.

    Parameters
    ----------
    G : nx.MultiDiGraph
        street graph

    Returns
    -------
    None

    """

    H = G.copy()

    for uvk, data in H.edges.items():

        u, v, k = uvk
        intermediary_nodes = data.get('_intermediary_nodes')
        if intermediary_nodes is None or len(intermediary_nodes) == 0:
            continue

        nodes = [u] + intermediary_nodes + [v]
        nodes_i = range(0, len(nodes))
        line = data['geometry']

        distances = np.linspace(0, 1, len(nodes))
        points = [line.interpolate(distance, normalized=True) for distance in distances]

        for i in nodes_i[1:-1]:
            point = points[i]
            G.add_node(nodes[i], x=point.x, y=point.y)

        for i in nodes_i[0:-1]:
            this_data = copy.copy(data)
            this_data['geometry'] = shapely.LineString((points[i], points[i+1]))
            this_data['length'] = None
            G.add_edge(nodes[i], nodes[i+1], **this_data)

        G.remove_edge(u, v, k)


def reset_intermediate_nodes(G):
    for uvk, data in G.edges.items():
        data['_intermediary_nodes'] = []
