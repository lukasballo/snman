import geopandas as gpd
import pandas as pd
import shapely as shp
import osmnx as ox
from . import io, geometry_tools
import networkx as nx
import itertools as it

def merge_nodes_geometric(G, tolerance, given_intersections_gdf=None, regions=None):
    """
    Geometrically merge nodes within some distance of each other.
    Parameters
    ----------
    G : networkx.MultiDiGraph
        a projected graph
    tolerance : float
        buffer nodes to this distance (in graph's geometry's units) then merge
        overlapping polygons into a single polygon via a unary union operation
    Returns
    -------
    merged : GeoSeries
        the merged overlapping polygons of the buffered nodes
    """
    if regions is None:
        # buffer nodes GeoSeries then get unary union to merge overlaps
        auto_intersections = ox.utils_graph.graph_to_gdfs(G, edges=False)["geometry"].buffer(tolerance).unary_union
    else:
        pass
        """
        auto_intersections = []
        for i, row in regions.iterrows():
            nodes = io._get_nodes_within_polygon(G, row['geometry'])
            H = G.subgraph(nodes)
            if len(H.nodes) == 0: continue
            sub_intersections = ox.utils_graph.graph_to_gdfs(H, edges=False)["geometry"].buffer(tolerance).unary_union
            sub_intersections = sub_intersections.intersection(row['geometry'])
            sub_intersections = geometry_tools.ensure_multipolygon(sub_intersections)
            sub_intersections = list(sub_intersections.geoms)
            intersections.extend(sub_intersections)
        print(intersections, len(intersections))
        intersections = geometry_tools.ensure_multipolygon(intersections)
        """

    if given_intersections_gdf is not None:
        given_intersections = geometry_tools.ensure_multipolygon(given_intersections_gdf['geometry'].unary_union)
        auto_intersections = geometry_tools.ensure_multipolygon(auto_intersections.difference(given_intersections))
        auto_intersections = gpd.GeoSeries(auto_intersections.geoms, crs=G.graph["crs"])
        auto_intersections_gdf = gpd.GeoDataFrame({'geometry': auto_intersections, 'point_geometry': auto_intersections.centroid})

    return pd.concat([given_intersections_gdf, auto_intersections_gdf])

def consolidate_intersections(G, intersections, reconnect_edges=True):
    """
    Consolidate intersections comprising clusters of nearby nodes.
    Merge nodes and return a rebuilt graph with consolidated intersections and
    reconnected edge geometries.
    The tolerance argument should be adjusted to approximately match street
    design standards in the specific street network, and you should always use
    a projected graph to work in meaningful and consistent units like meters.
    Returned graph's node IDs represent clusters rather than osmids. Refer to
    nodes' osmid_original attributes for original osmids. If multiple nodes
    were merged together, the osmid_original attribute is a list of merged
    nodes' osmids.
    Parameters
    ----------
    G : networkx.MultiDiGraph
        a projected graph
    tolerance : float
        nodes are buffered to this distance (in graph's geometry's units) and
        subsequent overlaps are dissolved into a single node
    reconnect_edges : bool
        ignored if rebuild_graph is not True. if True, reconnect edges and
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
    node_clusters = intersections
    centroids = node_clusters.centroid
    node_clusters["x"] = centroids.x
    node_clusters["y"] = centroids.y


    # STEP 2
    # attach each node to its cluster of merged nodes. first get the original
    # graph's node points then spatial join to give each node the label of
    # cluster it's within
    node_points = ox.utils_graph.graph_to_gdfs(G, edges=False)[["geometry", "street_count"]]
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
            # identify all the (weakly connected) component in cluster
            wccs = list(nx.weakly_connected_components(G.subgraph(nodes_subset.index)))
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
        if len(osmids) == 1:
            # if cluster is a single node, add that node to new graph
            osmid = osmids[0]
            H.add_node(cluster_label, osmid_original=osmid, **G.nodes[osmid])

        else:
            # if cluster is multiple merged nodes, create one new node to
            # represent them
            H.add_node(
                cluster_label,
                osmid_original=str(osmids),
                x=nodes_subset["x"].iloc[0],
                y=nodes_subset["y"].iloc[0],
            )

    # calculate street_count attribute for all nodes lacking it
    null_nodes = [n for n, sc in H.nodes(data="street_count") if sc is None]
    street_count = ox.stats.count_streets_per_node(H, nodes=null_nodes)
    nx.set_node_attributes(H, street_count, name="street_count")

    if not G.edges or not reconnect_edges:
        # if reconnect_edges is False or there are no edges in original graph
        # (after dead-end removed), then skip edges and return new graph as-is
        return H

    # STEP 6
    # create new edge from cluster to cluster for each edge in original graph
    gdf_edges = ox.utils_graph.graph_to_gdfs(G, nodes=False)
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