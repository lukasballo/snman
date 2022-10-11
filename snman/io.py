import networkx as nx
from . import osmnx as ox
from . import config, graph_tools
import geopandas as gpd
import pyproj
import shapely.ops
import momepy

def export_streetgraph_to_shp(street_graph):
    """
    Exports street graph as a shape file

    Params
    ------
    street_graph : nx.MultiGraph
    """
    nodes, edges = ox.graph_to_gdfs(street_graph)
    edges.columns = [column[0:10] for column in edges.columns]
    #edges['ln_desc'] = edges['ln_desc'].apply(lambda ln_desc: print(type(ln_desc)))
    edges['ln_desc'] = edges['ln_desc'].apply(lambda ln_desc: ' | '.join(ln_desc))
    """
    for ln_desc in edges['ln_desc']:
        print(ln_desc)
        ' | '.join(ln_desc)
    """
    #edges['ln_desc'] = edges['ln_desc'].apply(lambda ln_desc: '*')
    export_gdf_to_shp(edges, config.data_path + 'edges.shp')


def export_gdf_to_shp(gdf, file_path):
    gdf.to_file(file_path, schema=None)


def import_shp_to_gdf(file_path):
    edges = gpd.read_file(file_path)
    edges = edges.to_crs(2056)
    return edges


def convert_crs_of_street_graph(street_graph, to_crs):
    from_crs = pyproj.CRS(street_graph.graph['crs'])
    to_crs = pyproj.CRS('EPSG:' + to_crs)
    project = pyproj.Transformer.from_crs(from_crs, to_crs, always_xy=True).transform

    # Update the street_graph's metadata
    street_graph.graph["crs"] = 'EPSG:' + str(config.crs)

    # Transform the geometry of all edges
    for edge in street_graph.edges(data=True, keys=True):
        if "geometry" in edge[3]:
            edge[3]["geometry"] = shapely.ops.transform(project, edge[3]["geometry"])
