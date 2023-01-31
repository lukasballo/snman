# import osmnx as ox
# import matplotlib.pyplot as plt
# import fiona
import geopandas as gpd
import copy
import shapely as shp
import networkx as nx
import snman.osmnx_customized as oxc
import snman

print('Starting...')

# Constants
INTERSECTION_TOLERANCE = 10
inputs_path = 'C:/DATA/CLOUD STORAGE/polybox/Research/SNMan/SNMan Shared/inputs/'
export_path = 'C:/DATA/CLOUD STORAGE/polybox/Research/SNMan/SNMan Shared/qgis_previews/'
oxc.utils.config(useful_tags_way=snman.constants.OSM_TAGS)

# =====================================================================================
# LOAD DATA
# =====================================================================================

print('Load perimeters')
perimeters = snman.load_perimeters(inputs_path + 'perimeters/perimeters.shp')

print('Get data from OSM server')
G = oxc.graph_from_polygon(
    perimeters.loc['hardbruecke']['geometry'],
    custom_filter=snman.constants.OSM_FILTER,
    simplify=True,
    simplify_strict=False,
    retain_all=True,
    one_edge_per_direction=False,
)