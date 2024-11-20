from . import street_graph, lane_graph, space_allocation
from .constants import *

"""
class StreetSection(list):
    def __init__(self, lanes):
        self.extend(lanes)

    def set_from_string(self, street_section_string):
        pass
"""

StreetCrossSection = space_allocation.SpaceAllocation
LaneCrossSection = space_allocation.Lane
