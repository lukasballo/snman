import geopandas as gpd
import pandas as pd
import shapely as shp
from . import osmnx_customized as oxc
from . import io, geometry_tools
import networkx as nx
import itertools as it


def merge_nodes_geometric(G, tolerance, given_intersections_gdf=None, regions_gdf=None):
    """
    Create intersection geometries.
    Optionally use given intersection geometries to locally override the auto-detected results

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
                .buffer(list(map(lambda n_streets: 1 if n_streets<3 else tolerance,G_gdf.street_count)))
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
    G : nx.MultiDiGraph
        street graph
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

    # STEP 1
    # buffer nodes to passed-in distance and merge overlaps. turn merged nodes
    # into gdf and get centroids of each cluster as x, y
    node_clusters = intersections_gdf
    centroids = node_clusters.centroid
    node_clusters["x"] = centroids.x
    node_clusters["y"] = centroids.y


    # STEP 2
    # attach each node to its cluster of merged nodes. first get the original
    # graph's node points then spatial join to give each node the label of
    # cluster it's within
    node_points = oxc.utils_graph.graph_to_gdfs(G, edges=False)[["geometry", "street_count", "highway"]]
    gdf = gpd.sjoin(node_points, node_clusters, how="left", predicate="within")
    gdf = gdf.drop(columns="geometry").rename(columns={"index_right": "cluster"})

    # STEP 3
    # if a cluster contains multiple components (i.e., it's not connected)
    # move each component to its own cluster (otherwise you will connect
    # nodes together that are not truly connected, e.g., nearby deadends or
    # surface streets with bridge).

    groups = gdf.groupby("cluster")
    for cluster_label, nodes_subset in groups:
        if len(nodes_subset) > 1:
            Gs = G.subgraph(nodes_subset.index).copy()
            # Ignore pedestrian links for the detection of weakly connected components
            #for id, edge in Gs.copy().edges.items():
            #    if edge.get('highway') in {'path', 'footway', 'steps'}:
            #        Gs.remove_edge(*id)
            # identify all the (weakly connected) component in cluster
            wccs = list(nx.weakly_connected_components(Gs))
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
    gdf["cluster"] = gdf["cluster"].factorize()[0]

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
        traffic_signals = 1 * ('traffic_signals' in highway_tags)
        if len(osmids) == 1:
            # if cluster is a single node, add that node to new graph
            osmid = osmids[0]
            H.add_node(cluster_label,
                osmid_original=osmid,
                traffic_signals = traffic_signals,
                **G.nodes[osmid]
            )

        else:
            # if cluster is multiple merged nodes, create one new node to
            # represent them
            H.add_node(
                cluster_label,
                osmid_original=str(osmids),
                highway=str(highway_tags),
                traffic_signals=traffic_signals,
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
            H.add_edge(u2, v2, **data)

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
                old_coords = list(H.edges[u, v, k]["geometry"].coords)
                new_coords = xy + old_coords if cluster_label == u else old_coords + xy
                new_geom = shp.ops.LineString(new_coords)
                H.edges[u, v, k]["geometry"] = new_geom

                # update the edge length attribute, given the new geometry
                H.edges[u, v, k]["length"] = new_geom.length

    return H
