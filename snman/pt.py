import geopandas as gpd
import networkx as nx
import warnings


def match_pt(G, pt_network):
    """
    Matches the public transport network onto the street network.

    Parameters
    ----------
    G: nx.MultiGraph
        street graph
    pt_network: gpd.GeoDataFrame
        Linestrings of public transport routes, based on the open dataset of ZVV
        "Linien des öffentlichen Verkehrs (OGD) (kantonaler Datensatz)"
        https://data.stadt-zuerich.ch/dataset/ktzh_linien_des_oeffentlichen_verkehrs__ogd_,
        with following additional columns:
            * TYPE (Tram, Bus, Nightbus, Microbus)
            * ALIGNMENT (tunnel, <None>)
    """

    # Add buffers around the pt routes
    pt_network_buffers = pt_network.copy()
    pt_network_buffers.geometry = pt_network_buffers.geometry.buffer(15, resolution=16)

    for edge in G.edges(data=True, keys=True):

        edge_geometry = edge[3]['geometry']

        # Get all pt routes whose buffers intersect with this edge
        pt_routes = pt_network_buffers[
            pt_network_buffers.intersects(edge_geometry).tolist()
            and pt_network_buffers['ALIGNMENT'] != 'tunnel'
        ].copy()

        if len(pt_routes) == 0:
            continue

        # Calculate the overlapping length
        with warnings.catch_warnings():
            # we suppress warnings due to a known issue in the intersection function
            # see here: https://github.com/shapely/shapely/issues/1345
            warnings.simplefilter("ignore")
            pt_routes['intersection_length_prop'] = [
                edge_geometry.intersection(pt_route).length / edge_geometry.length if edge_geometry.length != 0 else 0
                for pt_route in pt_routes.geometry
            ]

        # Keep only those that overlap over a substantial part
        pt_routes = pt_routes[pt_routes['intersection_length_prop'] > 0.7]

        edge[3]['pt_routes'] = str(pt_routes.LINIENNUMM.tolist())
        # TODO: Distinguish directions of pt lines
        edge[3]['pt_tram'] = pt_routes['TYPE'].eq('Tram').any() * 1
        edge[3]['pt_bus'] = pt_routes['TYPE'].eq('Bus').any() * 1
        edge[3]['pt_mcrbus'] = pt_routes['TYPE'].eq('Microbus').any() * 1
