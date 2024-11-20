# this is an experimental setup for a new, unified data structure that will create a permanent sync between
# street graph and lane graph

from . import street_graph, lane_graph, cross_section
from .constants import *

class Network:
    def __init__(self):
        self.referenced_attributes_from_lanes_to_streets = []
        # Street Graph
        self.G = street_graph.StreetGraph()
        # A list of Lane Graphs, one for each allocation
        self.Ls = {}


    def add_node(self, n, **data):
        self.G.add_node(n, **data)
        # add the same node to all lane graphs but without the attributes
        for L in self.Ls:
            L.add_node(n)


    def add_street(self, u, v, **data):
        k = self.G.add_edge(u, v, **data)
        if u > v:
            self.reverse_street(u, v, k)
        return k


    def reverse_street(self, u, v, k):
        pass


    def get_street(self, u, v, k):
        data = self.G[(u,v,k)]
        return (u,v,k), data


    def remove_street(self, u, v, k):
        # remove all lanes on that street first
        for allocation, L in self.Ls.items():
            lanes = self.get_lanes_of_street(u, v, k, allocation)
            for uvl, data in lanes:
                L.remove_edge(uvl)
        # then remove the street
        self.G.remove_edge(u, v, k)


    def add_allocation(self, key):
        self.Ls[key] = lane_graph.LaneGraph()
        for n, data in self.G.nodes.items():
            self.Ls[key].add_node(n, **data)


    def add_lane(self, lu, lv, k, allocation, cross_section, **data):
        # translate the lane nodes to street nodes which are always in ascendant order
        u, v = sorted((lu, lv))
        # determine a unique lane id within the street
        existing_lanes = self.get_lanes_of_street(u, v, k, allocation)
        ids = [data['id'] for uvl, data in existing_lanes]
        if len(ids) > 0:
            new_id = max(ids) + 1
        else:
            new_id = 0
        # generate an l-key for the lane
        lk = f"{str(u)}-{str(v)}-{str(k)}-{str(new_id)}"
        # add the new lane into the lane graph
        self.Ls[allocation].add_edge(
            lu, lv, lk,
            u=u, v=v, k=k, id=new_id,
            cross_section=cross_section,
            **data
        )
        return lk


    def get_lane(self, lu, lv, lk, allocation):
        data = self.Ls[allocation].edges[(lu, lv, lk)]
        return (lu, lv, lk), data


    def remove_lane(self, lu, lv, lk, allocation):
        self.Ls[allocation].remove_edge(lu, lv, lk)


    def get_lanes_of_street(self, u, v, k, allocation):
        L = self.Ls[allocation]
        M = L.subgraph([u, v])
        edges = M.edges.items()
        # filter by street key to remove all lanes from other parallel streets
        edges = [(uvl, data) for uvl, data in edges if data['k'] == k]
        return edges


    def get_cross_section_of_street(self, u, v, k, allocation):
        # return an object derived from the currently used "SpaceAllocation"
        pass


    def set_lanes_of_street(self, cross_section_object):
        pass


    def reference_variables_from_lanes_to_streets(self, attrs):
        self.referenced_attributes_from_lanes_to_streets = attrs
        for allocation, L in self.Ls.items():
            for attribute in self.referenced_attributes_from_lanes_to_streets:
                for luvk, data in L.edges.items():
                    lu, lv, lk = luvk
                    u = data['u']
                    v = data['v']
                    k = data['k']
                    data[attribute] = lambda: self.G.edges[(u, v, k)][attribute]

