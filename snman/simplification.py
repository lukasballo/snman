import geopandas as gpd
import pandas as pd
import shapely as shp
from . import osmnx_customized as oxc
from . import io, geometry_tools, graph, street_graph, merge_edges, street_graph_node
from .constants import *
import networkx as nx
import itertools as it
import copy


def simplify_edge_geometries(G, radius=DEFAULT_SIMPLIFICATION_RADIUS):
    for uvk, edge in G.edges.items():
        if edge.get('geometry') is not None and edge.get('_include_in_simplification', True):
            edge['geometry'] = edge['geometry'].simplify(radius, preserve_topology=False)


def create_intersection_geometries(G, tolerance=DEFAULT_INTERSECTION_TOLERANCE, given_intersections_gdf=None, regions_gdf=None):
    """
    Create intersection geometries.
    Optionally use given intersection geometries to locally override the auto-detected results.

    Derived from merge_nodes_geometric() in osmnx.

    Parameters
    ----------
    G : nx.MultiGraph or nx.MultiDiGraph
        street graph
    tolerance : int
        radius of circles around the nodes
    given_intersections_gdf : gpd.GeoDataFrame
        explicitly defined intersection geometries,
        created by .io.load_intersections()
    regions_gdf : gpd.GeoDataFrame
        simplification regions,
        created by .io.load_regions()

    Returns
    -------
    intersections_gdf : gpd.GeoDataFrame
        updated intersection geometries,
        including a mixture of the auto-detected and explicitly defined intersections
    """

    G_gdf = oxc.utils_graph.graph_to_gdfs(G, edges=False)

    # exclude nodes already covered by given intersections
    # this is important to get rid of residuals from the buffers after subtracting the given intersections
    if given_intersections_gdf is not None:
        G_gdf = G_gdf[~G_gdf.within(given_intersections_gdf.unary_union)]

    if regions_gdf is None:
        # buffer nodes GeoSeries then get unary union to merge overlaps
        auto_intersections = (
            G_gdf["geometry"]
                .buffer(
                    list(
                        G_gdf.apply(
                            lambda row:
                                1 if row['street_count'] < 3 and len(row.get('layers', [])) == 1 else tolerance,
                            axis=1
                        )
                    )
                )
                .unary_union
        )
    else:
        # for every region, create buffers and clip them with the region polygon
        auto_intersections = regions_gdf.apply(
            lambda region: shp.geometry.MultiPolygon(
                G_gdf["geometry"]
                    # use a small buffer for those nodes that have less than 3 streets
                    # this way, we avoid the creation of large intersections through chaining of unnecessary nodes
                    .buffer(
                        list(
                            map(
                                lambda n_streets: 1 if n_streets<3 else region.tolerance,
                                G_gdf.street_count
                            )
                        )
                    )
                    .unary_union
                    .intersection(region.geometry)
            ),
            axis=1
        )
        # create a list of lists with individual polygons
        auto_intersections = list(map(lambda multipolygon: multipolygon.geoms, auto_intersections))
        # flatten into a single list of polygons
        auto_intersections = list(it.chain.from_iterable(auto_intersections))
        # convert into a single multipolygon (across all regions)
        auto_intersections = shp.geometry.MultiPolygon(auto_intersections)


    if given_intersections_gdf is not None:
        given_intersections = geometry_tools.ensure_multipolygon(given_intersections_gdf['geometry'].unary_union)
        # subtract the given intersection from every detected intersection separately to avoid a unary union of
        # the resulting geometries
        auto_intersections = list(
            map(lambda geom:
                geometry_tools.ensure_multipolygon(geom.difference(given_intersections)),
                list(auto_intersections.geoms)
            )
        )

    auto_intersections = gpd.GeoSeries(auto_intersections, crs=G.graph["crs"])

    intersections_gdf = gpd.GeoDataFrame({
        'geometry': auto_intersections
    })

    if given_intersections_gdf is not None:
        intersections_gdf = pd.concat([given_intersections_gdf, intersections_gdf], ignore_index=True)

    intersections_gdf['point_geometry'] = intersections_gdf.centroid

    return intersections_gdf


def consolidate_intersections(G, intersections_gdf, reconnect_edges=True):
    """
    Merge nodes into larger intersections using intersection geometries.

    This function is a further development of osmnx.consolidate_intersections.

    Parameters
    ----------
    G : nx.MultiGraph
        centerline graph
    intersections_gdf : gpd.GeoSeries
        intersection geometries
    reconnect_edges : bool
        if True, reconnect edges and
        their geometries in rebuilt graph to the consolidated nodes and update
        edge length attributes; if False, returned graph has no edges (which
        is faster if you just need topologically consolidated intersection
        counts).
    Returns
    -------
    H : networkx.MultiDiGraph
        a rebuilt graph with consolidated intersections and reconnected
        edge geometries
    """

    G = copy.deepcopy(G)

    # STEP 1
    # prepare the provided intersection polygons
    node_clusters = intersections_gdf
    centroids = node_clusters.centroid
    node_clusters["x"] = centroids.x
    node_clusters["y"] = centroids.y
    # shift the cluster ids to avoid collision with node ids
    node_clusters['cluster'] = node_clusters.index + max(G.nodes)


    # STEP 2
    # attach each node to its cluster of merged nodes. first get the original
    # graph's node points then spatial join to give each node the label of
    # cluster it's within
    node_points = oxc.utils_graph.graph_to_gdfs(G, edges=False)[
        ["geometry", "street_count", "highway", "traffic_signals", "_include_in_simplification"]
    ]
    gdf = gpd.sjoin(node_points, node_clusters, how="left", predicate="within")
    gdf = gdf.drop(columns="geometry")

    gdf['osmid_temp'] = gdf.index
    # overwrite cluster ids with the normal ids in case of nodes that are excluded from simplification
    gdf['cluster'] = gdf.apply(lambda row: row['cluster'] if row['_include_in_simplification'] else row['osmid_temp'], axis=1)
    gdf.drop(columns='osmid_temp')

    gdf['intersection_index'] = gdf['index_right']

    # STEP 3
    # if a cluster contains multiple components (i.e., it's not connected)
    # move each component to its own cluster (otherwise you will connect
    # nodes together that are not truly connected, e.g., nearby deadends or
    # surface streets with bridge).

    groups = gdf.groupby("cluster")
    for cluster_label, nodes_subset in groups:
        if len(nodes_subset) > 1:
            Gsub = G.subgraph(nodes_subset.index).copy()
            # Ignore pedestrian links for the detection of weakly connected components
            #for id, edge in Gs.copy().edges.items():
            #    if edge.get('highway') in {'path', 'footway', 'steps'}:
            #        Gs.remove_edge(*id)
            # identify all the (weakly connected) component in cluster
            wccs = list(nx.weakly_connected_components(Gsub))
            if len(wccs) > 0:
                # if there are multiple components in this cluster
                suffix = 0
                for wcc in wccs:
                    # set subcluster xy to the centroid of just these nodes
                    idx = list(wcc)
                    subcluster_centroid = node_points.loc[idx].unary_union.centroid
                    gdf.loc[idx, "x"] = subcluster_centroid.x
                    gdf.loc[idx, "y"] = subcluster_centroid.y
                    # move to subcluster by appending suffix to cluster label
                    gdf.loc[idx, "cluster"] = f"{cluster_label}-{suffix}"
                    suffix += 1

    # give nodes unique integer IDs (subclusters with suffixes are strings)
    gdf['cluster'] = gdf['cluster'].factorize()[0]

    # STEP 4
    # create new empty graph and copy over misc graph data
    H = nx.MultiDiGraph()
    H.graph = G.graph

    # STEP 5
    # create a new node for each cluster of merged nodes
    # regroup now that we potentially have new cluster labels from step 3
    groups = gdf.groupby("cluster")
    for cluster_label, nodes_subset in groups:

        osmids = nodes_subset.index.to_list()
        highway_tags = set(nodes_subset['highway'].to_list())
        intersection_index = list(nodes_subset['intersection_index'])[0]

        traffic_signals = 1 * (1 in set(nodes_subset.get('traffic_signals').to_list()))

        if len(osmids) == 1:
            # if cluster is a single node, add that node to new graph
            osmid = osmids[0]
            H.add_node(
                cluster_label,
                osmid_original=osmid,
                traffic_signals=traffic_signals,
                highway=G.nodes[osmid].get('highway'),
                _include_in_simplification=G.nodes[osmid].get('_include_in_simplification'),
                x=G.nodes[osmid].get('x'),
                y=G.nodes[osmid].get('y')
            )

        else:
            # if cluster is multiple merged nodes, create one new node to
            # represent them
            H.add_node(
                cluster_label,
                osmid_original=str(osmids),
                traffic_signals=traffic_signals,
                highway=str(highway_tags),
                _include_in_simplification=node_clusters['simplify'][intersection_index],
                x=nodes_subset["x"].iloc[0],
                y=nodes_subset["y"].iloc[0],
            )

    # calculate street_count attribute for all nodes lacking it
    null_nodes = [n for n, sc in H.nodes(data="street_count") if sc is None]
    street_count = oxc.stats.count_streets_per_node(H, nodes=null_nodes)
    nx.set_node_attributes(H, street_count, name="street_count")

    if not G.edges or not reconnect_edges:
        # if reconnect_edges is False or there are no edges in original graph
        # (after dead-end removed), then skip edges and return new graph as-is
        return H

    # STEP 6
    # create new edge from cluster to cluster for each edge in original graph
    gdf_edges = oxc.utils_graph.graph_to_gdfs(G, nodes=False)
    for u, v, k, data in G.edges(keys=True, data=True):
        u2 = gdf.loc[u, "cluster"]
        v2 = gdf.loc[v, "cluster"]

        # only create the edge if we're not connecting the cluster
        # to itself, but always add original self-loops
        if (u2 != v2) or (u == v):
            data["u_original"] = u
            data["v_original"] = v
            if "geometry" not in data:
                data["geometry"] = gdf_edges.loc[(u, v, k), "geometry"]
            key2 = H.add_edge(u2, v2, **data)

            # reverse the edge attributes if its topological direction has changed
            #if u2 > v2:
            #    street_graph.reverse_edge(H, u2, v2, key2, reverse_topology=True)

    # STEP 7
    # for every group of merged nodes with more than 1 node in it, extend the
    # edge geometries to reach the new node point
    for cluster_label, nodes_subset in groups:

        # but only if there were multiple nodes merged together,
        # otherwise it's the same old edge as in original graph
        if len(nodes_subset) > 1:

            # get coords of merged nodes point centroid to prepend or
            # append to the old edge geom's coords
            x = H.nodes[cluster_label]["x"]
            y = H.nodes[cluster_label]["y"]
            xy = [(x, y)]

            # for each edge incident on this new merged node, update its
            # geometry to extend to/from the new node's point coords
            in_edges = set(H.in_edges(cluster_label, keys=True))
            out_edges = set(H.out_edges(cluster_label, keys=True))
            for u, v, k in in_edges | out_edges:
                geometry = H.edges[u, v, k]["geometry"]

                # for longer geometries, strip a part of the edge geometry before connecting it to the new node
                if geometry.length > 40:
                    strip_length = 10
                    if cluster_label == u:
                        geometry = shp.ops.substring(geometry, strip_length, geometry.length)
                        new_coords = xy + list(geometry.coords)
                    else:
                        geometry = shp.ops.substring(geometry, 0, geometry.length - strip_length)
                        new_coords = list(geometry.coords) + xy

                # for shorter geometries, just take the straight line
                else:
                    old_coords = list(geometry.coords)
                    new_coords = xy + [old_coords[-1]] if cluster_label == u else [old_coords[0]] + xy

                new_geom = shp.ops.LineString(new_coords)
                H.edges[u, v, k]["geometry"] = new_geom
                H.edges[u, v, k]["_extension"] = str([u, v, k])

                # update the edge length attribute, given the new geometry
                H.edges[u, v, k]["length"] = new_geom.length

    return street_graph.StreetGraph(H)


def split_through_edges_in_intersections(G, intersections_gdf):
    """
    Within each intersection polygon, split edges that are passing through it without having a node there.
    This is helpful for a proper simplification of complex intersections

    Parameters
    ----------
    G : nx.MultiGraph
        centerline graph
    intersections_gdf : gpd.GeoDataFrame
        intersection geometries

    Returns
    -------
    None
    """

    # abbreviations used:
    # ix: intersection
    #  e: edge

    intersections_gdf = copy.deepcopy(intersections_gdf)
    # retain the geometry and centroid in separate attributes for later joining
    intersections_gdf['ix_geometry'] = intersections_gdf['geometry']

    edges = oxc.graph_to_gdfs(G, nodes=False)
    edges = edges[edges['_include_in_simplification'] == True]
    # retain the edge geometry in a separate attribute for later joining
    edges['e_geometry'] = edges['geometry']

    # create new column with multipoint holding endpoints of each edge,
    # we will need them to detect if an edge start/ends in an intersection buffer
    edges['e_endpoints'] = edges.apply(
        lambda edge: shp.ops.MultiPoint([
            edge['geometry'].coords[0],
            edge['geometry'].coords[-1]
        ])
        # make sure the edge geometry is not corrupted
        if isinstance(edge.get('geometry'), shp.ops.LineString)
            and edge.get('geometry').is_valid
            and not edge.get('geometry').is_empty
        else None,
        axis=1
    )

    # build pairs of intersections and intersecting edges, keep index from edges
    a = gpd.sjoin(edges, intersections_gdf, how="inner", predicate="intersects", lsuffix='e', rsuffix='i')
    # stop here if there are no edge/intersection pairs
    if len(a) == 0:
        return

    # add a new column telling us whether the edge starts/ends within the intersection
    a['split_point'] = a.apply(
        lambda x: x['ix_geometry'].intersection(x['e_geometry']).centroid,
        axis=1
    )
    # add a new column with a meaningful split point
    a['edge_endpoint_in_intersection'] = a.apply(
        lambda x: x['ix_geometry'].intersects(x['e_endpoints']),
        axis=1
    )

    # filter for only those intersection/edge pairs where the edge does not start/end in the intersection
    a = a[a['edge_endpoint_in_intersection'] == False]

    a = a.groupby(['u', 'v', 'key'])['split_point'].unique()
    a = pd.DataFrame({'split_points': a})

    # stop here if there are no instances to process
    if len(a) == 0:
        return

    a.apply(
        lambda row: street_graph.split_edge(G, *row.name, row.split_points),
        axis=1
    )


def connect_components_in_intersections(G, intersections_gdf, separate_layers=True):
    """
    Creates connections between weakly connected components so that they can be merged into a single intersection
    In this process, we use the following restrictions:
    - ignoring dead-end nodes (so that we don't connect adjacent dead-ends)
    - connecting only components that are on the same layer, e.g. on the same level of a multi-level intersection

    Parameters
    ----------
    G : nx.MultiGraph
        centerline graph
    intersections_gdf : gpd.GeoDataFrame
        intersection geometries
    separate_layers : bool
        avoid connecting components that are on different layers (e.g., intersections above each other)

    Returns
    -------
    None
    """

    # get the nodes as a geodataframe
    node_points = oxc.graph_to_gdfs(G, edges=False)[["geometry", "street_count", "highway", '_include_in_simplification']]
    # eliminate dead ends from the process to keep them as they are
    node_points = node_points.query('street_count != 1')
    # assign each node to an intersection (polygon)
    gdf = gpd.sjoin(node_points, intersections_gdf, how="left", predicate="within")
    # clean up the columns of the resulting geodataframe (cluster=id of the intersection geometry)
    gdf = gdf.drop(columns="geometry").rename(columns={"index_right": "cluster"})
    # shift the cluster ids to avoid collision with node ids
    gdf['cluster'] = gdf['cluster'] + max(G.nodes)

    gdf['osmid_temp'] = gdf.index
    # overwrite cluster ids with the normal ids in case of nodes that are excluded from simplification
    gdf['cluster'] = gdf.apply(lambda row: row['cluster'] if row['_include_in_simplification'] else row['osmid_temp'], axis=1)
    gdf.drop(columns='osmid_temp')

    groups = gdf.groupby("cluster")
    for cluster_label, nodes_subset in groups:
        if len(nodes_subset) > 1:
            # identify all (weakly connected) components in cluster
            # wccs is a list of sets containing id's of the nodes in each weakly connected component
            wccs = list(nx.weakly_connected_components(G.subgraph(nodes_subset.index)))
            # skip this node if there are not at least two components (one component = no need for further connections)
            if len(wccs) <= 1:
                continue
            nodes = []
            layers = {}
            # iterate over all combinations of the selected nodes (one from each wcc) to create a complete from-to mesh
            for combination in it.combinations(range(len(wccs)), 2):

                a = combination[0]
                b = combination[1]

                # create gdfs of nodes from each wcc
                a = oxc.graph_to_gdfs(G.subgraph(wccs[a]), edges=False)
                b = oxc.graph_to_gdfs(G.subgraph(wccs[b]), edges=False)

                # extract the layers of both wcc
                a_layers = set().union(*a['layers'].tolist())
                b_layers = set().union(*b['layers'].tolist())

                # join with the nearest points
                joined = gpd.sjoin_nearest(a, b, distance_col="distance")
                joined['node_a'] = joined.index
                joined = joined.rename(columns={"index_right": "node_b"})

                # get the closest pair
                closest_pair = joined.sort_values(by='distance', ascending=True).to_dict('records')[0]
                a = closest_pair['node_a']
                b = closest_pair['node_b']

                a_data = G.nodes[a]
                b_data = G.nodes[b]

                # skip this connector if the layer sets don't match
                if separate_layers and a_layers.isdisjoint(b_layers):
                    continue

                geom = shp.ops.LineString((
                    shp.ops.Point(a_data.get('x'), a_data.get('y')),
                    shp.ops.Point(b_data.get('x'), b_data.get('y'))
                ))

                G.add_edge(
                    a, b, geometry=geom,
                    osmid=0,
                    _components_connector=True
                )


def add_layers_to_nodes(G):
    """
    Add a set to each node representing the layers to which is belongs.

    Parameters
    ----------
    G : nx.MultiGraph or nx.MultiDiGraph
        street graph

    Returns
    -------
    None
    """
    for i, data in G.nodes.items():
        edges = list(G.in_edges(i, keys=True, data=True)) + list(G.out_edges(i, keys=True, data=True))
        layers = set([edge[3].get('layer', 0) for edge in edges])
        data['layers'] = layers


def simplify_street_graph(
        G,
        given_intersections_gdf,
        exclude_edge_hierarchies_from_simplification=[],
        edge_geometries_simplification_radius=None,
        iterations=4,
        intersection_tolerance=10,
        verbose=False
):

    #print('Light simplification of edge geometries')
    simplify_edge_geometries(G, radius=2)

    for i in range(iterations):
        if verbose:
            print('ITERATION', i)

        # here, you can choose whether some roads should be excluded from the simplification
        def add_include_in_simplification_to_edge(G, uvk):
            data = G.edges[uvk]
            data['_include_in_simplification'] = (
                    data.get('hierarchy') not in exclude_edge_hierarchies_from_simplification
            )

        graph.apply_function_to_each_edge(G, add_include_in_simplification_to_edge)

        def add_include_in_simplification_to_node(G, n):
            uvks = set(G.in_edges(n, keys=True)).union(set(G.out_edges(n, keys=True)))
            edges_included_in_simplification = [G.edges[uvk].get('_include_in_simplification') for uvk in uvks]
            # include node in simplification only if all of its edges are also included in simplification
            G.nodes[n]['_include_in_simplification'] = (
                    G.nodes[n].get('_include_in_simplification', True)
                    and False not in edges_included_in_simplification
            )

        graph.apply_function_to_each_node(G, add_include_in_simplification_to_node)

        #print('Detect intersections')
        add_layers_to_nodes(G)
        intersections_gdf = create_intersection_geometries(
            G,
            tolerance=intersection_tolerance,
            given_intersections_gdf=given_intersections_gdf
        )

        #print('Split through edges in intersections')
        split_through_edges_in_intersections(G, intersections_gdf)

        #print('Detect intersections (repeat to ensure that no points are outside of intersections)')
        add_layers_to_nodes(G)
        intersections_gdf = create_intersection_geometries(
            G,
            tolerance=intersection_tolerance,
            given_intersections_gdf=given_intersections_gdf
        )

        #print('Add layers and hierarchies to nodes')
        add_layers_to_nodes(G)
        graph.apply_function_to_each_node(G, street_graph_node.add_hierarchies)

        #print('Add connections between components in intersections')
        connect_components_in_intersections(G, intersections_gdf, separate_layers=True)

        #print('Consolidate intersections')
        G = consolidate_intersections(
            G, intersections_gdf,
            reconnect_edges=True
        )

        #print('Merge consecutive edges')
        add_layers_to_nodes(G)
        merge_edges.merge_consecutive_edges(G)

        #print('Merge parallel edges')
        merge_edges.merge_parallel_edges(G)

        #print('Update precalculated attributes')
        street_graph.update_precalculated_attributes(G)

        # In rare cases, the geometries get corrupted throughout the process.
        # In this step, we fill them with straight lines.
        street_graph.fill_wrong_edge_geometries(G)

    #print('Heavy simplification of edge geometries')
    if edge_geometries_simplification_radius:
        simplify_edge_geometries(G, radius=edge_geometries_simplification_radius)

    #print('Keep only the largest connected component')
    G = graph.keep_only_the_largest_connected_component(G, weak=True)

    return G, intersections_gdf