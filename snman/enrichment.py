from leuvenmapmatching.matcher.distance import DistanceMatcher
from leuvenmapmatching.map.inmem import InMemMap
import shapely as shp
from statistics import mean

def match_linestrings(G, source, source_to_target_columns):

    map_con = InMemMap("source", use_latlon=False, use_rtree=True, index_edges=True, crs_xy=2056)

    # please note that lv works with lat, lon (reverse order)
    for id, data in G.nodes.items():
        map_con.add_node(id, (data['y'], data['x']))

    for id, data in G.edges.items():
        u = int(id[0])
        v = int(id[1])
        # the graph is directed so we need to add an edge in each direction
        map_con.add_edge(u,v)
        map_con.add_edge(v,u)

    matcher = DistanceMatcher(map_con, max_dist=20, max_dist_init=20, max_lattice_width=5, non_emitting_states=True, only_edges=True)

    def _get_nodes_of_linestring(geom):
        if isinstance(geom, shp.geometry.MultiLineString):
            geom = geom.geoms[0]
        path = geom.coords
        path = [coords[::-1] for coords in path]
        matcher.match(path)
        nodes = matcher.path_pred_onlynodes
        #print(nodes)
        return nodes

    source['nodes'] = source.apply(lambda x: _get_nodes_of_linestring(x['geometry']), axis=1)

    for source_column, target_column in source_to_target_columns.items():

        edge_values = {}
        for index, edge in source.iterrows():
            value = edge[source_column]
            nodes = edge['nodes']
            if len(nodes) < 2:
                continue
            node_pairs = [nodes[i:i+2] for i in range(len(nodes)-1)]
            for node_pair in node_pairs:
                u = node_pair[0]
                v = node_pair[1]
                if edge_values.get((u,v)) is None:
                    edge_values[(u,v)] = []
                edge_values[(u,v)].append(value)

        for id, data in G.edges.items():
            data[target_column + '_forward'] = mean(edge_values.get(id[:2], [0]))
            data[target_column + '_backward'] = mean(edge_values.get(id[:2][::-1], [0]))