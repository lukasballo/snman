from .constants import *
from . import space_allocation, street_graph
from . import osmnx_customized as oxc
import copy


def street_sections(G, key_lanes_description=KEY_LANES_DESCRIPTION, weight='length'):

    G = copy.deepcopy(G)

    space_allocation.reorder_lanes(G, lanes_attribute=KEY_LANES_DESCRIPTION_AFTER)
    street_graph.organize_edge_directions(
        G,
        method='by_top_order_lanes',
        key_lanes_description=KEY_LANES_DESCRIPTION_AFTER
    )

    nodes, edges = oxc.graph_to_gdfs(G, nodes=True, edges=True)
    edges[key_lanes_description] = edges[key_lanes_description].apply(lambda x: '|'.join(x))
    result = edges[[key_lanes_description, weight]].groupby(key_lanes_description).sum(weight)
    result = result.sort_values(weight, ascending=False)

    return result


def street_sections_change(
        G,
        key_lanes_description=KEY_LANES_DESCRIPTION,
        key_lanes_description_after=KEY_LANES_DESCRIPTION_AFTER,
        weight='length'
):

    G = copy.deepcopy(G)



