import networkx as nx
from . import utils
from .constants import *


def cast_attributes_for_export(G, uvk, lane_keys):
    """
    Cast edge attributes to strings for export.

    Parameters
    ----------
    G : nx.MultiDiGraph
        Street graph
    uvk : tuple
        Edge identifier (u, v, key)
    lane_keys : list
        List of lane description keys to convert to strings
    """
    data = G.edges[uvk]

    data['sensors_forward'] = utils.safe_dumps(data.get('sensors_forward'))
    data['sensors_backward'] = utils.safe_dumps(data.get('sensors_backward'))
    data['_intermediary_nodes'] = utils.safe_dumps(data.get('_intermediary_nodes'))

    for lane_key in lane_keys:
        data[lane_key] = ' | '.join(utils.convert_list_items_to_strings(data.get(lane_key, [])))
