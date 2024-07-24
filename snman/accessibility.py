from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
import geopandas as gpd
from .constants import *
import datetime
import r5py
import shapely as shp
import math
import osmnx_customized as oxc
import networkx as nx


# run the mode choice model to calculate each person's generalized cost to each statent destination
def calculate_behavioral_cost(
        euclidean_distance,
        travel_time_private_cars=np.inf, travel_time_access_egress_private_cars=np.inf, path_length_private_cars=np.inf,
        travel_time_transit=np.inf,
        travel_time_cycling=np.inf, travel_time_access_egress_cycling=np.inf,
        travel_time_foot=np.inf, travel_time_access_egress_foot=np.inf,
        age=None
):
    """
    Calculates individual generalized cost based on the personal properties and mode choice.
    Based on mode choice model implemented in eqasim.

    Returns a tuple of mode specific choice probabilities (dict) and the resulting behavioral travel time (float)
    """

    U = {}

    X_dist = euclidean_distance / 1000
    Mi_hhIncome = 12260
    X_hhIncome = 10000
    X_work = 0.5
    X_city_center = 0.2
    Mi_dist = 39
    Lambda_distTT = 0.1147
    Beta_cost = -0.088
    Lambda_dist_cost = -0.2209
    Lambda_hhIncome = -0.8169

    X_IVT_car = travel_time_private_cars / 60
    X_AET_car = travel_time_access_egress_private_cars / 60
    X_cost_car = 0.26 * path_length_private_cars / 1000
    Alpha_car = -0.8
    Beta_TT_car = -0.0192
    Beta_TT_walk = -0.0457
    Beta_work_car = -1.1606
    Beta_city_center_car = -0.4590

    U[MODE_PRIVATE_CARS] = (
            Alpha_car
            + Beta_TT_car * X_IVT_car * (X_dist / Mi_dist) ** Lambda_distTT
            + Beta_TT_walk * X_AET_car
            + Beta_cost * (X_dist / Mi_dist) ** Lambda_dist_cost + X_cost_car * (
                        X_hhIncome / Mi_hhIncome) ** Lambda_hhIncome
            + Beta_work_car * X_work
            + Beta_city_center_car * X_city_center
    )

    X_railTT = travel_time_transit * 0.8 / 60  # surrogate
    X_busTT = travel_time_transit * 0.2 / 60  # surrogate
    X_AET_PT = 0  # = 0 because R5 includes it in the travel time
    X_waiting_time = 5  # surrogate
    X_number_of_connections = travel_time_transit / 60 / 20
    X_cost_PT = min(2.7, euclidean_distance / 1000 * 0.6)
    X_headway = 10  # surrogate
    Alpha_PT = 0
    Beta_railTT = -0.0072
    Beta_busTT = -0.0124
    Beta_AET_PT = -0.0142
    Beta_wait = -0.0124
    Beta_lineSwitch = -0.17
    Beta_headway = -0.0301
    Beta_OVGK = -0.8  # surrogate

    U[MODE_TRANSIT] = (
            Alpha_PT
            + (Beta_railTT * X_railTT + Beta_busTT * X_busTT) * (X_dist / Mi_dist) ** Lambda_distTT
            + Beta_AET_PT * X_AET_PT
            + Beta_wait * X_waiting_time
            + Beta_lineSwitch * X_number_of_connections
            + Beta_cost * (X_dist / Mi_dist) ** Lambda_dist_cost * X_cost_PT * (
                        X_hhIncome / Mi_hhIncome) ** Lambda_hhIncome
            + Beta_headway * X_headway
            + Beta_OVGK
    )

    X_bikeTT = travel_time_cycling / 60 + travel_time_access_egress_cycling / 60
    X_age = age
    Alpha_bike = -0.1258
    Beta_TT_bike = -0.1258
    Beta_age_60_plus_bike = -2.6588

    U[MODE_CYCLING] = (
            Alpha_bike
            + Beta_TT_bike * X_bikeTT + (X_dist / Mi_dist) ** Lambda_distTT
            + Beta_age_60_plus_bike * (X_age >= 60)
    )

    X_walkTT = travel_time_foot / 60 + travel_time_access_egress_foot / 60
    Alpha_walk = 0.5903
    Beta_TT_walk = -0.0457
    Theta_threshold_walkTT = 120

    U[MODE_FOOT] = (
            Alpha_walk
            + Beta_TT_walk * X_walkTT * (X_dist / Mi_dist) ** Lambda_distTT
            + (1 - 100 ** (X_walkTT / Theta_threshold_walkTT))
    )

    # replace nan with -inf for non-existent options (=infinite cost)
    U = {mode: -np.inf if np.isnan(U_mode) else U_mode for mode, U_mode in U.items()}

    # mode-specific choice probabilities
    denominator = sum([np.exp(U_mode) for U_mode in U.values()])
    if denominator != 0:

        P = {mode: np.exp(U_mode) / denominator for mode, U_mode in U.items()}
        # print(U)

        # mode-specific cost
        C = {mode: -U_mode for mode, U_mode in U.items()}

        # weighted cost
        C_weighted = sum([C[mode] * P[mode] if not np.isinf(C[mode]) else 0 for mode in P.keys()])

    else:

        P = {mode: None for mode, U_mode in U.items()}
        C = {mode: np.inf for mode, U_mode in U.items()}
        C_weighted = np.inf

    return P, C, C_weighted

















def calculate_accessibility_for_statent_cell(
        r5_transit_network,
        G_modes,
        L_modes,
        statpop_with_statent_ids,
        statent,
        statent_id,
        departure_time,
        distance_limit=10000,
        destinations_sample=1,
        population_sample=1,
        sum_accessibility_contributions=True,
        cost_function_exponent=-0.7,
        cycling_speed_kmh=18,
        walking_speed_kmh=1.34*3.6,
        access_egress_detour_factor=1.5,
        min_travel_time=5*60,
        min_euclidian_distance=100
):
    """
    Calculates accessibility for every resident associated with a given cell in the statent dataset.
    Using R5 for transit and networkx for other modes.

    Parameters
    ----------
    r5_transit_network: r5py.TransportNetwork
        transport network without travel time adjustments
    G_modes: dict
        pre-filtered street graphs for cycling, foot, and private cars
    L_modes: dict
        pre-calculated street graphs for cycling, foot, and private cars
    statpop_with_statent_ids: gpd.GeoDataFrame
        statpop dataset, with the id's of closest statent cells
    statent: gpd.GeoDataFrame
        statent dataset
    statent_id: int
        which cell should be calculated, we can't calculate all cells at once due to memory constraints and because
        r5 does not allow a cutoff distance
    distance_limit: int
        cutoff crowfly distance
    min_travel_time: float
        minimum travel time in seconds to avoid 0 travel times
    min_euclidian_distance: float
        minimum Euclidean distance in meters to avoid 0 distances


    Returns
    -------
    gpd.GeoDataFrame
    """

    # --------------------------------------------

    # filter statpop for residents within the given origin statent cell
    statpop_filtered = statpop_with_statent_ids.query(f'statent_id == {statent_id}')
    statpop_filtered = statpop_filtered.sample(frac=population_sample, random_state=0)
    statpop_filtered

    # --------------------------------------------

    if len(statpop_filtered) == 0:
        return

    # --------------------------------------------

    # make a light version of statpop
    statpop_reduced = statpop_filtered.reset_index()[[
        'record', 'age', 'sex', 'maritalstatus', 'residencepermit', 'residentpermit', 'statent_id', 'geometry'
    ]]
    statpop_reduced.set_index('record', inplace=True)
    statpop_reduced

    # --------------------------------------------

    # filter statent for cells within a distance limit from the origin
    origin = statent.set_index('id').geometry[statent_id]
    origins = gpd.GeoDataFrame(
        {
            "id": [statent_id],
            "geometry": [origin]
        },
        crs="EPSG:2056",
    )
    statent_filtered = statent[statent['geometry'].within(
        origin.buffer(distance_limit)
    )]
    statent_filtered = gpd.GeoDataFrame(statent_filtered)
    statent_filtered

    # --------------------------------------------

    # make a sample of statpop cells but inlcude the origin cell in every case
    statent_filtered = pd.concat([
        statent_filtered[statent_filtered['id'] == statent_id],
        statent_filtered[statent_filtered['id'] != statent_id].sample(frac=destinations_sample, random_state=0)
    ])

    # --------------------------------------------

    # match statent points to lanegraph nodes
    # TODO: do this only once

    for mode in [MODE_PRIVATE_CARS, MODE_CYCLING, MODE_FOOT]:
        # match statent points to nearest nodes
        nodes = oxc.nearest_nodes(
            L_modes[mode],
            list(map(lambda geom: geom.x, statent_filtered.geometry)),
            list(map(lambda geom: geom.y, statent_filtered.geometry)),
            return_dist=True
        )

        statent_filtered[[f'closest_node_{mode}', f'egress_{mode}']] = list(zip(*nodes))

    statent_filtered

    # --------------------------------------------

    tt_matrices = {}

    # --------------------------------------------

    # calculate transit travel time matrix using R5
    tt_computer_transit = r5py.TravelTimeMatrixComputer(
        r5_transit_network,
        transport_modes=[r5py.TransportMode.TRANSIT],
        origins=origins,
        destinations=statent_filtered,
        departure=departure_time,
        snap_to_network=True
    )
    tt_matrices[MODE_TRANSIT] = tt_computer_transit.compute_travel_times()
    tt_matrices[MODE_TRANSIT].rename(columns={
        'travel_time': 'travel_time_transit',
        'from_id': 'from_cell',
        'to_id': 'to_cell'
    }, inplace=True)
    # convert from minutes to seconds
    tt_matrices[MODE_TRANSIT]['travel_time_transit'] *= 60
    tt_matrices[MODE_TRANSIT]['travel_time_egress_transit'] = 0

    tt_matrices[MODE_TRANSIT]

    # --------------------------------------------

    # calculate shortest paths to all destinations in networkx
    for mode in [MODE_PRIVATE_CARS, MODE_FOOT, MODE_CYCLING]:

        print(mode)

        origin_cell = statent_filtered[statent_filtered['id'] == statent_id]
        origin_node = min(origin_cell[f'closest_node_{mode}'])
        access_cost = min(origin_cell[f'egress_{mode}'])
        print(origin_node, access_cost)

        if mode == MODE_PRIVATE_CARS:
            weight = 'matsim_tt_cars'
        else:
            weight = f'cost_{mode}'

        # calculate cost to all other nodes
        dijkstra_path_costs = nx.single_source_dijkstra_path_length(L_modes[mode], origin_node, weight=weight)

        # convert the dijkstra path costs dict to a dataframe like from R5
        tt_matrices[mode] = pd.DataFrame(list(dijkstra_path_costs.items()), columns=['to_node', f'path_length_{mode}'])
        tt_matrices[mode]['from_node'] = origin_node

    tt_matrices[MODE_PRIVATE_CARS]

    # --------------------------------------------

    # merge the tt matrix of transit on cell_ids
    statent_filtered_2 = pd.merge(
        statent_filtered, tt_matrices[MODE_TRANSIT],
        left_on='id', right_on='to_cell', how='left'
    )
    statent_filtered_2.drop(columns=['from_cell', 'to_cell'], inplace=True)
    statent_filtered_2[f'travel_time_total_{MODE_TRANSIT}'] = statent_filtered_2[f'travel_time_{MODE_TRANSIT}']
    statent_filtered_2

    # --------------------------------------------

    # merge the other tt matrices on node ids
    for mode in [MODE_PRIVATE_CARS, MODE_FOOT, MODE_CYCLING]:

        # merge the travel time matrices with the statent dataset
        # such that every destination has the cost of reaching it
        statent_filtered_2 = pd.merge(
            statent_filtered_2, tt_matrices[mode],
            left_on=f'closest_node_{mode}', right_on='to_node', how='left'
        )
        statent_filtered_2.drop(columns=['from_node', 'to_node'], inplace=True)

        origin_cell = statent_filtered[statent_filtered['id'] == statent_id]
        access_cost = min(origin_cell[f'egress_{mode}'])

        if mode == MODE_FOOT:
            speed_factor = walking_speed_kmh / 3.6
        elif mode == MODE_CYCLING:
            speed_factor = cycling_speed_kmh / 3.6
        else:
            speed_factor = 1

        # calculate total travel time by adding access and egress cost and considering speed factors
        statent_filtered_2[f'travel_time_{mode}'] = (
            (statent_filtered_2[f'path_length_{mode}'] / speed_factor)
        )

        statent_filtered_2[f'travel_time_access_egress_{mode}'] = (
                + ((access_cost * access_egress_detour_factor) / (walking_speed_kmh / 3.6))
                + ((statent_filtered_2[f'egress_{mode}'] * access_egress_detour_factor) / (walking_speed_kmh / 3.6))
        )

    statent_filtered_2

    # --------------------------------------------

    statent_filtered_2['from_cell'] = statent_id
    statent_filtered_2

    # --------------------------------------------

    # add euclidean distance

    origin_cell = statent_filtered[statent_filtered['id'] == statent_id]
    origin_geometry = min(origin_cell.geometry)
    statent_filtered_2['euclidean_distance'] = statent_filtered_2.apply(
        lambda row: shp.distance(origin_geometry, row.geometry),
        axis=1
    )
    statent_filtered_2

    # --------------------------------------------

    destinations_with_cost = pd.merge(
        statpop_reduced.reset_index(),
        statent_filtered_2.reset_index()[[
            'id', 'VOLLZEITAEQ_TOTAL',
            'euclidean_distance',
            'travel_time_foot',
            'travel_time_cycling',
            'travel_time_transit',
            'travel_time_private_cars',
            'travel_time_access_egress_foot',
            'travel_time_access_egress_cycling',
            'travel_time_access_egress_private_cars',
            'path_length_private_cars',
            'path_length_cycling',
            'path_length_foot',
            'closest_node_private_cars',
            'from_cell',
            'geometry'
        ]].rename(columns={'geometry': 'destination_geometry'}),
        how='left', left_on='statent_id', right_on='from_cell')
    destinations_with_cost

    # --------------------------------------------

    # ensure minimum travel time and euclidian distance
    for mode in [MODE_CYCLING, MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_FOOT]:
        destinations_with_cost[f'travel_time_{mode}'] = np.maximum(destinations_with_cost[f'travel_time_{mode}'],
                                                                   min_travel_time)

    destinations_with_cost['euclidean_distance'] = np.maximum(destinations_with_cost['euclidean_distance'],
                                                              min_euclidian_distance)

    # --------------------------------------------

    destinations_with_cost[['choice_probabilities', 'cost', 'weighted_cost']] = destinations_with_cost.apply(
        lambda row: calculate_behavioral_cost(
            **row[[
                'euclidean_distance',
                'travel_time_private_cars', 'travel_time_access_egress_private_cars', 'path_length_private_cars',
                'travel_time_transit',
                'travel_time_cycling', 'travel_time_access_egress_cycling',
                'travel_time_foot', 'travel_time_access_egress_foot',
                'age'
            ]]
        ),
        axis=1,
        result_type='expand'
    )

    destinations_with_cost

    # --------------------------------------------

    destinations_with_cost = pd.concat([
        destinations_with_cost,
        pd.json_normalize(destinations_with_cost['choice_probabilities']).add_prefix('p_'),
        pd.json_normalize(destinations_with_cost['cost']).add_prefix('c_'),
    ],
        axis=1
    )

    # calculate total accessibility contribution using cost function and destination utility
    destinations_with_cost['accessibility_contribution'] = (
            (destinations_with_cost['weighted_cost'] ** -0.7)
            * destinations_with_cost['VOLLZEITAEQ_TOTAL']
    )

    destinations_with_cost

    # --------------------------------------------

    # calculate accessibility contribution for every mode separately
    for mode in [MODE_CYCLING, MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_FOOT]:
        destinations_with_cost[f'accessibility_contribution_{mode}'] = (
                (destinations_with_cost[f'c_{mode}'] ** -0.7)
                * destinations_with_cost['VOLLZEITAEQ_TOTAL']
        )

    destinations_with_cost

    # --------------------------------------------

    # for each person, sum the accessibility contributions across all destinations
    accessibility = destinations_with_cost.groupby('record').agg({
        'age': 'first',
        'sex': 'first',
        'maritalstatus': 'first',
        'residencepermit': 'first',
        'residentpermit': 'first',
        'statent_id': 'first',
        'accessibility_contribution': 'sum',
        'accessibility_contribution_cycling': 'sum',
        'accessibility_contribution_foot': 'sum',
        'accessibility_contribution_private_cars': 'sum',
        'accessibility_contribution_transit': 'sum',
        'geometry': 'first'
    })

    accessibility.rename(columns={
        'accessibility_contribution': 'accessibility',
        'accessibility_contribution_cycling': 'accessibility_cycling',
        'accessibility_contribution_foot': 'accessibility_foot',
        'accessibility_contribution_private_cars': 'accessibility_private_cars',
        'accessibility_contribution_transit': 'accessibility_transit',
    }, inplace=True)

    accessibility.reset_index(inplace=True)
    accessibility = gpd.GeoDataFrame(accessibility, crs=2056)
    accessibility

    # --------------------------------------------

    return accessibility
