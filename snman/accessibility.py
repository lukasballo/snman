from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
import geopandas as gpd
from .constants import *
import datetime
import r5py
import shapely as shp
import math
from snman import osmnx_customized as oxc
import networkx as nx


# run the mode choice model to calculate each person's generalized cost to each statent destination
def calculate_behavioral_cost(
        euclidean_distance,
        travel_time_private_cars=np.inf, travel_time_access_egress_private_cars=np.inf, path_length_private_cars=np.inf,
        travel_time_transit=np.inf,
        travel_time_cycling=np.inf, travel_time_access_egress_cycling=np.inf, path_length_cycling=np.inf,
        travel_time_pedelec=np.inf, travel_time_access_egress_pedelec=np.inf, path_length_pedelec=np.inf,
        travel_time_s_pedelec=np.inf, travel_time_access_egress_s_pedelec=np.inf, path_length_s_pedelec=np.inf,
        travel_time_foot=np.inf, travel_time_access_egress_foot=np.inf,
        sex=None,
        age=None
):
    """
    old mode choice model
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
    X_cost_car = 0.26 * euclidean_distance / 1000 * 1.2  # surrogate
    Alpha_car = -0.8
    Beta_TT_car = -0.0192
    Beta_TT_walk = -0.0457
    Beta_work_car = -1.1606
    Beta_city_center_car = -0.4590

    U[MODE_PRIVATE_CARS] = (
            Alpha_car
            + Beta_TT_car * X_IVT_car * (X_dist / Mi_dist) ** Lambda_distTT
            + Beta_TT_walk * X_AET_car
            + Beta_cost
                * (X_dist / Mi_dist) ** Lambda_dist_cost
                * X_cost_car * (X_hhIncome / Mi_hhIncome) ** Lambda_hhIncome
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
            + Beta_cost
                * (X_dist / Mi_dist) ** Lambda_dist_cost
                * X_cost_PT * (X_hhIncome / Mi_hhIncome) ** Lambda_hhIncome
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

    U[MODE_PEDELEC] = -np.inf
    U[MODE_S_PEDELEC] = -np.inf

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

    C_logsum_exponent = np.sum([np.exp(v_ijm) for v_ijm in U.values()])

    return P, C, C_weighted, C_logsum_exponent








# with new mode choice model
def calculate_behavioral_cost_updated(
        euclidean_distance,
        travel_time_private_cars=np.inf, travel_time_access_egress_private_cars=np.inf, path_length_private_cars=np.inf,
        travel_time_transit=np.inf,
        travel_time_cycling=np.inf, travel_time_access_egress_cycling=np.inf, path_length_cycling=np.inf,
        travel_time_pedelec=np.inf, travel_time_access_egress_pedelec=np.inf, path_length_pedelec=np.inf,
        travel_time_s_pedelec=np.inf, travel_time_access_egress_s_pedelec=np.inf, path_length_s_pedelec=np.inf,
        travel_time_foot=np.inf, travel_time_access_egress_foot=np.inf, path_length_foot=np.inf,
        sex=None,
        age=None,
):
    """
    Calculates individual generalized cost based on the personal properties and mode choice.
    Based on 2024-12 "E-bike City" Mode choice model, Meyer de Freitas et al.

    Returns a tuple of mode specific choice probabilities (dict) and the resulting behavioral travel time (float)
    """

    U = {}

    Alpha_car = 0.3
    Beta_TT_car = -6.0
    Beta_parking_cost = -0.305164
    Beta_cost = -0.06934
    Beta_externalities_car = 0.644314
    Gamma_externalitiesByKm_car = 0.1601
    X_in_vehicle_distance_car = euclidean_distance * 1.2 / 1000  # surrogate, we don't know the real distance
    X_cost_car = 0.188 * euclidean_distance * 1.2 / 1000
    X_TT_car = travel_time_private_cars / 3600 + travel_time_access_egress_private_cars / 3600
    X_parking_cost = 4  # surrogate

    U[MODE_PRIVATE_CARS] = (
            Alpha_car
            + Beta_TT_car * X_TT_car
            + Beta_parking_cost * X_parking_cost
            + Beta_cost * X_cost_car
            + Beta_externalities_car * Gamma_externalitiesByKm_car * X_in_vehicle_distance_car
    )


    Alpha_PT = -0.8
    Beta_female_PT = 0.31345
    Beta_age_PT = 0.00354
    Beta_degurba2_PT = -0.94476
    Beta_degurba3_PT = -1.25242
    Beta_TT_PT = -2.0
    Beta_access_egress_time_PT = -1.96973
    Beta_freq = -0.50346
    Beta_cost = -0.06934
    Beta_externalities_PT = 1.44709
    Gamma_externalitiesByKm_PT = 0.08
    X_sex_female = sex==2
    X_age = age
    X_degurba_medium = 0.33             # surrogate, 1/3
    X_degurba_low = 0.33                # surrogate, 1/3
    X_TT_PT = travel_time_transit / 3600
    X_access_egress_time_PT = 0         # R5 includes it in the travel time
    X_freq_PT = 6                       # surrogate, 10 min headway
    X_PT_dist = 1.5 * euclidean_distance / 1000   # surrogate, assuming 1.5 detour factor
    X_in_vehicle_distance_PT = X_PT_dist
    X_cost_PT = 0.5 * 0.5 * max([
        3.4,
        (X_PT_dist <= 5) * 0.89 * X_PT_dist
        + (X_PT_dist >= 5) * 0.589 * X_PT_dist
    ])

    U[MODE_TRANSIT] = (
            Alpha_PT
            + Beta_female_PT * X_sex_female
            + Beta_age_PT * X_age
            + Beta_degurba2_PT * X_degurba_medium
            + Beta_degurba3_PT * X_degurba_low
            + Beta_TT_PT * X_TT_PT
            + Beta_access_egress_time_PT * X_access_egress_time_PT
            + Beta_freq * X_freq_PT
            + Beta_cost * X_cost_PT
            + Beta_externalities_PT * Gamma_externalitiesByKm_PT * X_in_vehicle_distance_PT
    )

    Alpha_bike = 1.4
    Beta_female_bike = -0.03147
    Beta_age_bike = -0.02074
    Beta_degurba2_bike = -1.29194
    Beta_degurba3_bike = -1.92303
    Beta_TT_bike = -2.4
    Beta_cost = -0.06934
    Beta_externalities_bike = 3.18593
    Gamma_externalitiesByKm_bike = -0.0364
    X_TT_bike = travel_time_cycling / 3600 + travel_time_access_egress_cycling / 3600
    X_cost_bike = 0
    X_distance_bike = path_length_cycling / 1000

    U[MODE_CYCLING] = (
            Alpha_bike
            + Beta_female_bike * X_sex_female
            + Beta_age_bike * X_age
            + Beta_degurba2_bike * X_degurba_medium
            + Beta_degurba3_bike * X_degurba_low
            + Beta_TT_bike * X_TT_bike
            + Beta_cost * X_cost_bike
            + Beta_externalities_bike * Gamma_externalitiesByKm_bike * X_distance_bike
    )

    Alpha_ebike = -0.8
    Beta_female_ebike = 0.32921
    Beta_age_ebike = 0.00268
    Beta_degurba2_ebike = -0.51416
    Beta_degurba3_ebike = -0.64266
    Beta_TT_ebike = -2
    Beta_cost = -0.06934
    Beta_externalities_ebike = 0
    Gamma_externalitiesByKm_ebike = 0.0264
    X_cost_ebike = 0
    X_TT_ebike = travel_time_pedelec / 3600  + travel_time_access_egress_pedelec / 3600

    U[MODE_PEDELEC] = (
            Alpha_ebike
            + Beta_female_ebike * X_sex_female
            + Beta_age_ebike * X_age
            + Beta_degurba2_ebike * X_degurba_medium
            + Beta_degurba3_ebike * X_degurba_low
            + Beta_TT_ebike * X_TT_ebike
            + Beta_cost * X_cost_ebike
            #+ Beta_externalities_ebike * Gamma_externalitiesByKm_ebike * X_distance_ebike    # no externalities
    )

    Alpha_spedelec = -1.4
    Beta_female_spedelec = -0.363751
    Beta_age_spedelec = -0.026566
    Beta_degurba2_spedelec = -0.464518
    Beta_degurba3_spedelec = -0.651147
    Beta_TT_spedelec = -0.3
    Beta_cost = -0.06934
    Beta_externalities_ebike = 0
    Gamma_externalitiesByKm_ebike = 0.0264
    X_cost_spedelec = 0
    X_TT_spedelec = travel_time_s_pedelec / 3600 + travel_time_access_egress_s_pedelec / 3600

    U[MODE_S_PEDELEC] = (
            Alpha_spedelec
            + Beta_female_spedelec * X_sex_female
            + Beta_age_spedelec * X_age
            + Beta_degurba2_spedelec * X_degurba_medium
            + Beta_degurba3_spedelec * X_degurba_low
            + Beta_TT_spedelec * X_TT_spedelec
            + Beta_cost * X_cost_spedelec
            #+ Beta_externalities_ebike * Gamma_externalitiesByKm_ebike * X_distance_ebike    # no externalities
    )

    Alpha_walk = 0.8
    Beta_female_walk = 0.07087
    Beta_age_walk = -0.00447
    Beta_degurba2_walk = -0.70167
    Beta_degurba3_walk = -0.37095
    Beta_TT_walk = -2
    Beta_cost = -0.06934
    Beta_externalities_walk = 0
    Gamma_externalitiesByKm_walk = -0.0997
    X_TT_walk = travel_time_foot / 3600  + travel_time_access_egress_foot / 3600
    X_cost_walk = 0

    U[MODE_FOOT] = (
            Alpha_walk
            + Beta_female_walk * X_sex_female
            + Beta_age_walk * X_age
            + Beta_degurba2_walk * X_degurba_medium
            + Beta_degurba3_walk * X_degurba_low
            + Beta_TT_walk * X_TT_walk
            + Beta_cost * X_cost_walk
            #+ Beta_externalities_walk * Gamma_externalitiesByKm_walk * X_distance_walk      # no externalities
    )

    # cut off cycling length at 40 km
    if path_length_cycling > 40 * 1000:
        U[MODE_CYCLING] = -np.inf

    # cut off cycling length at 40 km
    if path_length_pedelec > 40 * 1000:
        U[MODE_PEDELEC] = -np.inf

    # cut off cycling length at 40 km
    if path_length_s_pedelec > 40 * 1000:
        U[MODE_S_PEDELEC] = -np.inf

    # cut off walking length at 5 km
    if path_length_foot > 5 * 1000:
        U[MODE_FOOT] = -np.inf

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

    C_logsum_exponent = np.sum([np.exp(v_ijm) for v_ijm in U.values()])

    return P, C, C_weighted, C_logsum_exponent
















def calculate_accessibility_for_statent_cell_logsum(
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
        cycling_speed_kmh=15,
        pedelec_speed_kmh=20,
        s_pedelec_speed_kmh=24,
        walking_speed_kmh=1.34*3.6,
        access_egress_detour_factor=1.5,
        min_travel_time=1*60,
        min_euclidian_distance=100,
        accessibility_beta = {
            MODE_PRIVATE_CARS: 1,
            MODE_TRANSIT: 1.4,
            MODE_CYCLING: 1.6,
            MODE_PEDELEC: 1.6,
            MODE_S_PEDELEC: 1.6,
            MODE_FOOT: 2
        },
        return_destinations_with_cost=False,
        mode_choice_model='ebc',
        cumulative_accessibility_cutoff_cost = 30 * 60
):
    """
    Calculates accessibility for every resident associated with a given cell in the statent dataset.
    Using R5 for transit and networkx for other modes.

    An updated, "proper" logsum definition.

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

    # make a sample of statent cells but include the origin cell in every case
    statent_filtered = pd.concat([
        statent_filtered[statent_filtered['id'] == statent_id],
        statent_filtered[statent_filtered['id'] != statent_id].sample(
            frac=destinations_sample,
            #weights='VOLLZEITAEQ_TOTAL',
            random_state=0
        )
    ])

    # --------------------------------------------

    # match statent points to lanegraph nodes
    # TODO: do this only once

    for mode in [MODE_PRIVATE_CARS, MODE_CYCLING, MODE_FOOT, MODE_PEDELEC, MODE_S_PEDELEC]:
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
    for mode in [MODE_PRIVATE_CARS, MODE_FOOT, MODE_CYCLING, MODE_PEDELEC, MODE_S_PEDELEC]:

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
    for mode in [MODE_PRIVATE_CARS, MODE_FOOT, MODE_CYCLING, MODE_PEDELEC, MODE_S_PEDELEC]:

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
        elif mode == MODE_PEDELEC:
            speed_factor = pedelec_speed_kmh / 3.6
        elif mode == MODE_S_PEDELEC:
            speed_factor = s_pedelec_speed_kmh / 3.6
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

        statent_filtered_2[f'travel_time_access_egress_transit'] = 0

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
            'travel_time_pedelec',
            'travel_time_s_pedelec',
            'travel_time_transit',
            'travel_time_private_cars',
            'travel_time_access_egress_foot',
            'travel_time_access_egress_cycling',
            'travel_time_access_egress_pedelec',
            'travel_time_access_egress_s_pedelec',
            'travel_time_access_egress_private_cars',
            'travel_time_access_egress_transit',
            'path_length_private_cars',
            'path_length_cycling',
            'path_length_pedelec',
            'path_length_s_pedelec',
            'path_length_foot',
            'closest_node_private_cars',
            'from_cell',
            'geometry'
        ]].rename(columns={'geometry': 'destination_geometry'}),
        how='left', left_on='statent_id', right_on='from_cell')
    destinations_with_cost

    # --------------------------------------------

    # ensure minimum travel time and euclidian distance
    for mode in [MODE_PRIVATE_CARS, MODE_CYCLING, MODE_FOOT, MODE_PEDELEC, MODE_S_PEDELEC]:
        destinations_with_cost[f'travel_time_{mode}'] = np.maximum(destinations_with_cost[f'travel_time_{mode}'],
                                                                   min_travel_time)

    destinations_with_cost['euclidean_distance'] = np.maximum(destinations_with_cost['euclidean_distance'],
                                                              min_euclidian_distance)

    # --------------------------------------------

    if mode_choice_model == 'ebc':

        destinations_with_cost[['choice_probabilities', 'cost', 'weighted_cost', 'logsum_exponent']] = destinations_with_cost.apply(
            lambda row: calculate_behavioral_cost_updated(
                **row[[
                    'euclidean_distance',
                    'travel_time_private_cars', 'travel_time_access_egress_private_cars', 'path_length_private_cars',
                    'travel_time_transit',
                    'travel_time_cycling', 'travel_time_access_egress_cycling', 'path_length_cycling',
                    'travel_time_pedelec', 'travel_time_access_egress_pedelec', 'path_length_pedelec',
                    'travel_time_s_pedelec', 'travel_time_access_egress_s_pedelec', 'path_length_s_pedelec',
                    'travel_time_foot', 'travel_time_access_egress_foot', 'path_length_foot',
                    'sex',
                    'age'
                ]],
            ),
            axis=1,
            result_type='expand'
        )

    else:

        destinations_with_cost[['choice_probabilities', 'cost', 'weighted_cost', 'logsum_exponent']] = destinations_with_cost.apply(
            lambda row: calculate_behavioral_cost(
                **row[[
                    'euclidean_distance',
                    'travel_time_private_cars', 'travel_time_access_egress_private_cars', 'path_length_private_cars',
                    'travel_time_transit',
                    'travel_time_cycling', 'travel_time_access_egress_cycling', 'path_length_cycling',
                    'travel_time_pedelec', 'travel_time_access_egress_pedelec', 'path_length_pedelec',
                    'travel_time_s_pedelec', 'travel_time_access_egress_s_pedelec', 'path_length_s_pedelec',
                    'travel_time_foot', 'travel_time_access_egress_foot', 'path_length_foot',
                    'sex',
                    'age'
                ]],
            ),
            axis=1,
            result_type='expand'
        )

    # Here we get the utilities * -1
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
    #destinations_with_cost['accessibility_contribution'] = (
    #        (destinations_with_cost['weighted_cost'] ** accessibility_beta)
    #        * destinations_with_cost['VOLLZEITAEQ_TOTAL']
    #)

    destinations_with_cost

    # --------------------------------------------

    # calculate accessibility contribution for every mode separately
    for mode in [MODE_CYCLING, MODE_PEDELEC, MODE_S_PEDELEC, MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_FOOT]:
        destinations_with_cost[f'accessibility_contribution_{mode}'] = (
                (destinations_with_cost[f'c_{mode}'] ** -accessibility_beta[mode])
                * (destinations_with_cost['VOLLZEITAEQ_TOTAL'] / destinations_sample)
        )

    # calculate cumulative accessibility contribution for every mode separately
    for mode in [MODE_CYCLING, MODE_PEDELEC, MODE_S_PEDELEC, MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_FOOT]:
        cost_col = f'c_{mode}'
        out_col = f'cumulative_accessibility_contribution_{mode}'

        mask = destinations_with_cost[cost_col] <= cumulative_accessibility_cutoff_cost
        destinations_with_cost[out_col] = 0

        # Only where the cost is below cutoff (i.e. within 30 minutes), use the number of jobs
        # as accessibility contribution. Everywhere else it will remain 0.
        destinations_with_cost.loc[mask, out_col] = destinations_with_cost.loc[mask, 'VOLLZEITAEQ_TOTAL']
        destinations_with_cost[out_col] = destinations_with_cost[out_col] / destinations_sample

    destinations_with_cost


    # scale the logsum exponent with destination opportunities
    destinations_with_cost['logsum_exponent'] = (
            destinations_with_cost['logsum_exponent']
            * (destinations_with_cost['VOLLZEITAEQ_TOTAL'] / destinations_sample)
    )

    destinations_with_cost

    # --------------------------------------------

    # scale the logsum exponent with destination opportunities for every mode separately
    for mode in [MODE_CYCLING, MODE_PEDELEC, MODE_S_PEDELEC, MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_FOOT]:
        destinations_with_cost[f'logsum_exponent_{mode}'] = (
                np.exp(-destinations_with_cost[f'c_{mode}'])
                * (destinations_with_cost['VOLLZEITAEQ_TOTAL'] / destinations_sample)
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
        'logsum_exponent': 'sum',
        'logsum_exponent_cycling': 'sum',
        'logsum_exponent_pedelec': 'sum',
        'logsum_exponent_s_pedelec': 'sum',
        'logsum_exponent_foot': 'sum',
        'logsum_exponent_private_cars': 'sum',
        'logsum_exponent_transit': 'sum',
        #'accessibility_contribution': 'sum',
        'accessibility_contribution_cycling': 'sum',
        'accessibility_contribution_pedelec': 'sum',
        'accessibility_contribution_s_pedelec': 'sum',
        'accessibility_contribution_foot': 'sum',
        'accessibility_contribution_private_cars': 'sum',
        'accessibility_contribution_transit': 'sum',
        'cumulative_accessibility_contribution_cycling': 'sum',
        'cumulative_accessibility_contribution_pedelec': 'sum',
        'cumulative_accessibility_contribution_s_pedelec': 'sum',
        'cumulative_accessibility_contribution_foot': 'sum',
        'cumulative_accessibility_contribution_private_cars': 'sum',
        'cumulative_accessibility_contribution_transit': 'sum',
        'geometry': 'first'
    })

    accessibility.rename(columns={
        'logsum_exponent': 'sum_of_logsum_exponents_all',
        'logsum_exponent_cycling': 'sum_of_logsum_exponents_cycling',
        'logsum_exponent_pedelec': 'sum_of_logsum_exponents_pedelec',
        'logsum_exponent_s_pedelec': 'sum_of_logsum_exponents_s_pedelec',
        'logsum_exponent_foot': 'sum_of_logsum_exponents_foot',
        'logsum_exponent_private_cars': 'sum_of_logsum_exponents_private_cars',
        'logsum_exponent_transit': 'sum_of_logsum_exponents_transit',
        #'accessibility_contribution': 'accessibility',
        'accessibility_contribution_cycling': 'hansen_accessibility_cycling',
        'accessibility_contribution_pedelec': 'hansen_accessibility_pedelec',
        'accessibility_contribution_s_pedelec': 'hansen_accessibility_s_pedelec',
        'accessibility_contribution_foot': 'hansen_accessibility_foot',
        'accessibility_contribution_private_cars': 'hansen_accessibility_private_cars',
        'accessibility_contribution_transit': 'hansen_accessibility_transit',
        'cumulative_accessibility_contribution_cycling': 'cumulative_accessibility_cycling',
        'cumulative_accessibility_contribution_pedelec': 'cumulative_accessibility_pedelec',
        'cumulative_accessibility_contribution_s_pedelec': 'cumulative_accessibility_s_pedelec',
        'cumulative_accessibility_contribution_foot': 'cumulative_accessibility_foot',
        'cumulative_accessibility_contribution_private_cars': 'cumulative_accessibility_private_cars',
        'cumulative_accessibility_contribution_transit': 'cumulative_accessibility_transit'
    }, inplace=True)

    accessibility['hansen_accessibility'] = (
        accessibility['hansen_accessibility_cycling']
        + accessibility['hansen_accessibility_pedelec']
        + accessibility['hansen_accessibility_s_pedelec']
        + accessibility['hansen_accessibility_foot']
        + accessibility['hansen_accessibility_private_cars']
        + accessibility['hansen_accessibility_transit']
    )

    for mode in ['all', MODE_CYCLING, MODE_PEDELEC, MODE_S_PEDELEC, MODE_FOOT, MODE_PRIVATE_CARS, MODE_TRANSIT]:
        accessibility[f'logsum_{mode}'] = np.log(accessibility[f'sum_of_logsum_exponents_{mode}'])

    accessibility.reset_index(inplace=True)
    accessibility = gpd.GeoDataFrame(accessibility, crs=2056)
    accessibility

    # --------------------------------------------

    if return_destinations_with_cost is True:
        return destinations_with_cost
    else:
        return accessibility
