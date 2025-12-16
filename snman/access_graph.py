import copy
import shapely
import networkx as nx
import math
from . import graph, utils
from .constants import *
import geopandas as gpd
from . import osmnx_customized as oxc


def get_edges_by_node(G, i):
    """
    Get all edges connected to a node.

    Parameters
    ----------
    G : nx.Graph
        Graph
    i : Any
        Node identifier

    Returns
    -------
    dict
        Dictionary mapping edge identifiers to edge data
    """
    uvks = G.edges(i)
    edges = {uvk: G.edges[uvk] for uvk in uvks}
    return edges


def get_present_parking_spots(L, uvk, vehicle_length):
    """
    Calculate number of parking spots that can fit along an edge.

    Parameters
    ----------
    L : nx.MultiDiGraph
        Lane graph
    uvk : tuple
        Edge identifier (u, v, key)
    vehicle_length : float
        Length of a single vehicle

    Returns
    -------
    float
        Number of parking spots
    """
    data = L.edges[uvk]
    length = data['length']
    return length / vehicle_length


def create_access_graph(
        L, access_needs, offstreet_parking=None,
        radius=100,
        lanetype=LANETYPE_PARKING_PARALLEL, vehicle_length=PARALLEL_PARKING_CAR_LENGTH
):
    """
    Create an access graph for parking spot assignment.

    Parameters
    ----------
    L : nx.MultiDiGraph
        Lane graph
    access_needs : gpd.GeoDataFrame
        GeoDataFrame with locations needing parking spots
    offstreet_parking : gpd.GeoDataFrame, optional
        GeoDataFrame with offstreet parking locations
    radius : float, optional
        Search radius for matching parking spots (default: 100)
    lanetype : str, optional
        Type of parking lane (default: LANETYPE_PARKING_PARALLEL)
    vehicle_length : float, optional
        Length of a single vehicle (default: PARALLEL_PARKING_CAR_LENGTH)

    Returns
    -------
    AccessGraph
        Access graph for parking assignment
    """
    #access_needs_buffered = copy.deepcopy(access_needs)
    #access_needs_buffered.geometry = access_needs.buffer(radius)
    #access_needs_buffered['need_id'] = access_needs_buffered.index
    A = AccessGraph(crs=L.graph['crs'], lanetype=lanetype, vehicle_length=vehicle_length)

    # add residential locations as nodes
    access_needs.apply(
        lambda x: A.add_node(
            x.name,
            geometry=x['geometry'],
            x=x['geometry'].x,
            y=x['geometry'].y,
            type='needs_parking_spots',
            parking_spots=x['parking_spots_needed']
        ),
        axis=1
    )

    # add offstreet parking as nodes
    if offstreet_parking is not None:
        offstreet_parking.apply(
            lambda x: A.add_node(
                f"offstreet_parking_{x.name}",
                geometry=x['geometry'],
                x=x['geometry'].x,
                y=x['geometry'].y,
                type=x['type'],
                parking_spots=x['parking_spots']
            ),
            axis=1
        )

    # add parking lanes as nodes
    for uvk, data in L.edges.items():
        if data['lanetype'] == lanetype:
            point = data['geometry'].interpolate(0.5, normalized=True)
            A.add_node(
                uvk,
                geometry=point,
                x=point.x,
                y=point.y,
                type='has_parking_spots',
                parking_spots=None
            )
            A.update_lane(L, uvk)

    # build possible location-street pairs
    #L_edges = oxc.utils_graph.graph_to_gdfs(L, nodes=False).query(f'lanetype=="{lanetype}"')

    A_nodes = oxc.utils_graph.graph_to_gdfs(A, edges=False)
    has_parking_spots = A_nodes.query("type=='has_parking_spots'")
    has_parking_spots['supply_id'] = has_parking_spots.index

    needs_parking_spots = A_nodes.query("type=='needs_parking_spots'")
    needs_parking_spots.geometry = needs_parking_spots.buffer(radius)
    needs_parking_spots['need_id'] = needs_parking_spots.index
    needs_parking_spots['parking_spots_needed'] = needs_parking_spots['parking_spots']

    joined = gpd.sjoin(needs_parking_spots.reset_index(), has_parking_spots.reset_index(), how='inner', predicate='intersects')
    #joined['uvk'] = joined.apply(
    #    lambda x: tuple(x[['u', 'v', 'key']]),
    #    axis=1
    #)
    pairs = joined[['supply_id', 'need_id', 'parking_spots_needed']]

    #print(pairs)

    # add location-street pairs
    pairs.reset_index().apply(
        lambda x: A.add_edge(
            x['need_id'], x['supply_id'],
            needs_parking_spots=x['need_id'],
            has_parking_spots=x['supply_id'],
            distance=shapely.distance(A.nodes[x['need_id']]['geometry'], A.nodes[x['supply_id']]['geometry']),
            parking_spots_assigned=x['parking_spots_needed']
        ),
        axis=1
    )

    return A


def gravity_model(A, iterations=15):
    """
    Apply gravity model to assign parking spots.

    Parameters
    ----------
    A : AccessGraph
        Access graph
    iterations : int, optional
        Number of iterations for the gravity model (default: 15)

    Returns
    -------
    None
    """
    # apply cost function on edges
    for uv, data in A.edges.items():
        data['f_cost'] = math.exp(-0.1 * data['distance'])

    # initialize origin factors (will move all under-/over-assignment to destinations)
    for i, data in A.nodes.items():
        if data['type'] == 'has_parking_spots':
            data['a_i'] = 1

    for iteration in range(iterations):

        # destination factors
        for j, j_data in A.nodes.items():
            if j_data['type'] == 'needs_parking_spots':
                a_j = 0
                edges = get_edges_by_node(A, j)
                for uv, edge_data in edges.items():
                    i_data = A.nodes[edge_data['has_parking_spots']]
                    a_j += i_data['a_i'] * i_data['parking_spots'] * edge_data['f_cost']
                    j_data['a_j'] = utils.safe_division(1, a_j)

        # origin factors
        for i, i_data in A.nodes.items():
            if i_data['type'] == 'has_parking_spots':
                a_i = 0
                edges = get_edges_by_node(A, i)
                for uv, edge_data in edges.items():
                    j_data = A.nodes[edge_data['needs_parking_spots']]
                    a_i += j_data['a_j'] * j_data['parking_spots'] * edge_data['f_cost']
                    i_data['a_i'] = utils.safe_division(1, a_i)

    for uv, edge_data in A.edges.items():
        i = edge_data['has_parking_spots']
        j = edge_data['needs_parking_spots']
        i_data = A.nodes[i]
        j_data = A.nodes[j]
        edge_data['parking_spots_assigned'] = i_data['a_i'] * i_data['parking_spots'] * j_data['a_j'] * j_data[
            'parking_spots'] * edge_data['f_cost']

    # Calculate under-/ over-assignment
    for i, data in A.nodes.items():
        data['_assigned_parking_spots'] = A.get_assigned_parking_spots(i)
        data['_parking_overassignment'] = data['_assigned_parking_spots'] - data['parking_spots']

    def sum_over_attribute_of_neighbors(A, i, attribute):
        neighbor_ids = A.neighbors(i)
        neighbor_overassignment = [A.nodes[neighbor_id]['_parking_overassignment'] for neighbor_id in neighbor_ids]
        return sum(neighbor_overassignment)

    for i, data in A.nodes.items():
        data['_overassignment_in_neighbors'] = sum_over_attribute_of_neighbors(A, i, '_parking_overassignment')


def effect_of_parking_removal_on_underassignment(A, ii, iterations=15):
    """
    Calculate the effect of removing parking spots on underassignment.

    Parameters
    ----------
    A : AccessGraph
        Access graph
    ii : list
        List of node identifiers to remove
    iterations : int, optional
        Number of iterations for gravity model (default: 15)

    Returns
    -------
    float
        Effect on underassignment (0 if effect is minimal)
    """
    B = copy.deepcopy(
        A.subgraph(set(A.nodes).difference(set(ii)))
    )

    gravity_model(B, iterations=iterations)

    underassignment_A = sum([min(data['_parking_overassignment'], 0) for uv, data in A.nodes.items()])
    underassignment_B = sum([min(data['_parking_overassignment'], 0) for uv, data in B.nodes.items()])

    effect = underassignment_B - underassignment_A
    if effect > -0.5:
        return 0
    else:
        return effect


class AccessGraph(nx.Graph, graph.SNManGenericGraph):

    def __init__(self, *args, lanetype=None, vehicle_length=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.lanetype = lanetype
        self.vehicle_length = vehicle_length

    def update_lane(self, L, uvk):
        """
        Update parking spots for a lane.

        Parameters
        ----------
        L : nx.MultiDiGraph
            Lane graph
        uvk : tuple
            Edge identifier
        """
        data = self.nodes[uvk]
        data['parking_spots'] = get_present_parking_spots(L, uvk, self.vehicle_length)

    def get_assigned_parking_spots(self, i):
        """
        Get total assigned parking spots for a node.

        Parameters
        ----------
        i : Any
            Node identifier

        Returns
        -------
        float
            Total assigned parking spots
        """
        pairs = get_edges_by_node(self, i)
        return sum([pair['parking_spots_assigned'] for pair in pairs.values()])

    def get_total_parking_spots(self):
        """
        Get total available parking spots in the graph.

        Returns
        -------
        float
            Total parking spots
        """
        return sum(
            [data['parking_spots'] for uv, data in self.nodes.items() if data['type'] == 'has_parking_spots']
        )

    def get_total_needed_parking_spots(self):
        """
        Get total needed parking spots in the graph.

        Returns
        -------
        float
            Total needed parking spots
        """
        return sum(
            [data['parking_spots'] for uv, data in self.nodes.items() if data['type'] == 'needs_parking_spots']
        )
