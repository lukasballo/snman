import networkx as nx
from . import utils
from .constants import *


def cast_attributes_for_export(G, uvk, lane_keys):

    data = G.edges[uvk]

    data['sensors_forward'] = utils.safe_dumps(data.get('sensors_backward'))
    data['sensors_backward'] = utils.safe_dumps(data.get('sensors_backward'))
    data['_intermediary_nodes'] = utils.safe_dumps(data.get('sensors_backward'))

    for lane_key in lane_keys:
        data[lane_key] = ' | '.join(utils.convert_list_items_to_strings(data.get(lane_key, [])))
