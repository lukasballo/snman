
"""
nodes = copy.deepcopy(street_graph.nodes(data=True))
for node in nodes:
    print(node)
    node_data = node[1]
    buffer = 5
    point = shp.ops.Point(node_data.get('x'), node_data.get('y'))

    edges = copy.deepcopy(street_graph.edges(data=True, keys=True))
    for edge in edges:
        edge_data = edge[3]
        line = edge_data.get('geometry')
        if line == None or line.distance(point) > buffer:
            continue
        if point.distance(shp.ops.Point(line.coords[0])) <= buffer \
                or point.distance(shp.ops.Point(line.coords[1])) <= buffer:
            continue
        nearest_point = shp.ops.nearest_points(point, line)[1]
        #print(edge)
        snman.graph_tools._split_edge(street_graph, edge, nearest_point)

"""