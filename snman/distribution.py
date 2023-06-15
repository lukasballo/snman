import networkx as nx
from . import constants, space_allocation, hierarchy
from .constants import *


def set_given_lanes(
        G,
        source_lanes_attribute=constants.KEY_LANES_DESCRIPTION,
        target_lanes_attribute=constants.KEY_GIVEN_LANES_DESCRIPTION,
        maintain_motorized_access_on_street_level = False,
        hierarchies_to_fix = set()
):
    """
    Sets which lanes are given due to external policy definitions
    e.g. dedicated lanes for public transport, bidirectional lanes for cars, etc.

    Parameters
    ----------
    G : nx.MultiGraph
        street graph

    Returns
    -------
    None
    """

    for id, data in G.edges.items():

        source_lanes = data[source_lanes_attribute]
        target_lanes = []

        lane_stats = space_allocation._lane_stats(source_lanes)

        # for street with hierarchy to fix keep everything as it is
        if data.get('hierarchy') in hierarchies_to_fix:
            target_lanes = source_lanes

        # for streets with public transit...
        elif data.get('pt_tram') or data.get('pt_bus'):
            if lane_stats.n_lanes_dedicated_pt_backward + lane_stats.n_lanes_dedicated_pt_both_ways == 0:
                # ...keep backward lane if there is no dedicated backward/both ways lane
                target_lanes += [LANETYPE_MOTORIZED + DIRECTION_BACKWARD]
            if lane_stats.n_lanes_dedicated_pt_both_ways + lane_stats.n_lanes_dedicated_pt_backward \
                    + lane_stats.n_lanes_dedicated_pt_forward > 0:
                if maintain_motorized_access_on_street_level:
                    target_lanes += [LANETYPE_MOTORIZED + DIRECTION_TBD]
                else:
                    target_lanes += [LANETYPE_MOTORIZED + DIRECTION_TBD_OPTIONAL]
            if lane_stats.n_lanes_dedicated_pt_forward + lane_stats.n_lanes_dedicated_pt_both_ways == 0:
                # ...keep forward lane if there is no dedicated forward/both ways lane
                target_lanes += [LANETYPE_MOTORIZED + DIRECTION_FORWARD]

        # for normal streets, keep one single-direction lane
        elif data.get('hierarchy') in [hierarchy.MAIN_ROAD, hierarchy.LOCAL_ROAD, hierarchy.HIGHWAY]:
            if maintain_motorized_access_on_street_level:
                target_lanes += [LANETYPE_MOTORIZED + DIRECTION_TBD]
            else:
                target_lanes += [LANETYPE_MOTORIZED + DIRECTION_TBD_OPTIONAL]

        data[target_lanes_attribute] = target_lanes
