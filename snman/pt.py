import math
import networkx as nx
from . import constants
import geopandas as gpd
import numpy as np
import networkx as nx


def match_pt(street_graph, pt_network):
    """
    Matches the public transport network onto the street network.

    Parameters
    ----------
    street_graph: nx.MultiGraph
        (u, v, key, data)
    pt_network: gpd.GeoDataFrame
        Linestrings of public transport routes, based on the open dataset of ZVV
        "Linien des öffentlichen Verkehrs (OGD) (kantonaler Datensatz)"
        https://data.stadt-zuerich.ch/dataset/ktzh_linien_des_oeffentlichen_verkehrs__ogd_,
        with following additional columns:
            - TYPE (Tram, Bus, Nightbus, Microbus)
            - ALIGNMENT (tunnel, <None>)
    """

    # Add buffers around the pt routes
    pt_network_buffers = pt_network.copy()
    pt_network_buffers.geometry = pt_network_buffers.geometry.buffer(15, resolution=16)



    for edge in street_graph.edges(data=True, keys=True):

        edge_geometry = edge[3]['geometry']

        # Get all pt routes whose buffers intersect with this edge
        pt_routes = pt_network_buffers[
            pt_network_buffers.intersects(edge_geometry).tolist()
            and pt_network_buffers['ALIGNMENT'] != 'tunnel'
        ].copy()

        if len(pt_routes) == 0:
            continue

        # Calculate the overlapping length
        pt_routes['intersection_length_prop'] = [
            edge_geometry.intersection(pt_route).length / edge_geometry.length
            for pt_route in pt_routes.geometry
        ]

        # Keep only those that overlap over a substantial part
        pt_routes = pt_routes[pt_routes['intersection_length_prop'] > 0.7]

        edge[3]['pt_routes'] = str(pt_routes.LINIENNUMM.tolist())
        # TODO: Distinguish directions of pt lines
        edge[3]['pt_tram'] = pt_routes['TYPE'].eq('Tram').any() * 1
        edge[3]['pt_bus'] = pt_routes['TYPE'].eq('Bus').any() * 1
        edge[3]['pt_mcrbus'] = pt_routes['TYPE'].eq('Microbus').any() * 1
