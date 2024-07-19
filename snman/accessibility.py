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


def calculate_behavioral_cost(private_cars=999999, transit=999999, cycling=999999, foot=999999, age=None):
    """
    Calculates individual generalized cost based on the personal properties and mode choice.
    Based on mode choice model in

    Hörl, S., M. Balać and K.W. Axhausen (2019)
    Pairing discrete mode choice models and agent-based transport simulation with MATSim,
    paper presented at the 98th Annual Meeting of the Transportation Research Board (TRB 2019),
    Washington DC, January 2019.

    Returns a tuple of mode specific choice probabilities (dict) and the resulting behavioral travel time (float)

    Parameters
    ----------
    private_cars: float
        travel time in minutes
    transit: float
        travel time in minutes
    cycling: float
        travel time in minutes
    foot: float
        travel time in minutes
    age: int
        age in years

    Returns
    -------
    (dict, float)
    """

    # mode-specific costs (travel times in minutes)
    tt = {
        MODE_PRIVATE_CARS: private_cars,
        MODE_TRANSIT: transit,
        MODE_CYCLING: cycling,
        MODE_FOOT: foot
    }

    for mode, cost in tt.items():
        if np.isnan(cost) or np.isinf(cost):
            # tt[mode] = np.inf
            tt[mode] = 99999999

    # mode choice model parameters from Hörl et al. (2019), pp. 11, Table 1
    # t: time

    alpha_car = 0.827  # [-]
    beta_tt_car = -0.0667  # [min-1]

    alpha_pt = 0.0  # [-]
    beta_transfer = -0.17  # [-]
    beta_in_vehicle_t = -0.0192  # [min-1]
    beta_transfer_t = -0.0384  # [min-1]
    beta_access_egress_t = -0.0804  # [min-1]

    alpha_bike = -0.1  # [-]
    beta_tt_bike = -0.0805  # [min-1]
    beta_age_bike = -0.0496  # [years-1]

    alpha_walk = 0.63  # [-]
    beta_tt_walk = -0.141  # [min-1]

    beta_cost = -0.126  # [CHF-1]
    lambda_factor = -0.4  # [-1]
    phi_avg_distance = 40.0  # [km]

    phi_parking_search_penalty = 6.0  # [min]
    phi_access_egress_walk_t = 5.0  # [min]

    # surrogate parameters to replace missing information (R5 provides only travel time)
    avg_speed_car = 25.0  # [kmh]
    detour_parameter_car = 1.3  # [-]
    x_crowfly_distance = tt[MODE_PRIVATE_CARS] / 60 * avg_speed_car / detour_parameter_car
    x_cost_car = tt[MODE_PRIVATE_CARS] * 0.2
    x_cost_pt = tt[MODE_TRANSIT] * 0.2

    x_transfers = round(tt[MODE_TRANSIT] / 30) if tt[MODE_TRANSIT] is not np.inf else np.inf
    x_transfer_t = x_transfers * 5
    x_access_egress_t = 5

    # mode-specific utilities from Hörl et al. (2019), pp. 10, eq. 3-6
    u = {
        MODE_PRIVATE_CARS:
            alpha_car
            + beta_tt_car * tt[MODE_PRIVATE_CARS]
            + beta_tt_car * phi_parking_search_penalty
            + beta_tt_walk * phi_access_egress_walk_t
            + beta_cost * (x_crowfly_distance / phi_avg_distance) ** lambda_factor * x_cost_car
        ,
        MODE_TRANSIT:
            alpha_pt
            + beta_transfer * x_transfers
            + beta_transfer_t * x_transfer_t
            + beta_access_egress_t * x_access_egress_t
            + beta_cost * (x_crowfly_distance / phi_avg_distance) ** lambda_factor * x_cost_pt
        ,
        MODE_CYCLING:
            alpha_bike
            + beta_tt_bike * tt[MODE_CYCLING]
            + beta_age_bike * max(0, age - 18)
        ,
        MODE_FOOT:
            alpha_walk
            + beta_tt_walk * tt[MODE_FOOT]
    }

    # mode-specific choice probabilities
    denominator = sum([np.exp(u_mode) for u_mode in u.values()])
    P = {mode: np.exp(u_mode) / denominator for mode, u_mode in u.items()}

    # total cost, after weighting the mode-specific cost by mode-specific choice probabilities
    tt_total = sum([tt[mode] * P[mode] for mode in tt.keys()])
    return P, tt_total








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
        access_egress_detour_factor=2,
        base_travel_time=5*60
):
    """
    Calculates accessibility for every resident associated with a given cell in the statent dataset.
    Using R5 for transit and networkx for other modes.

    Parameters
    ----------
    transport_network_default: r5py.TransportNetwork
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
    base_travel_time: float
        travel time in minutes to be added to all travel timed, a small number is useful to prevent 0 min. travel times
    sum_accessibility_contributions: bool
        if True, the accessibility contributions will be summed into one accessibility value per person

    Returns
    -------
    gpd.GeoDataFrame
    """

    # filter statpop for residents within the given origin statent cell
    statpop_filtered = statpop_with_statent_ids.query(f'statent_id == {statent_id}')
    statpop_filtered = statpop_filtered.sample(frac=population_sample, random_state=0)
    statpop_filtered

    if len(statpop_filtered) == 0:
        return

    # make a light version of statpop
    statpop_reduced = statpop_filtered.reset_index()[[
        'record', 'age', 'sex', 'maritalstatus', 'residencepermit', 'residentpermit', 'statent_id', 'geometry'
    ]]
    statpop_reduced.set_index('record', inplace=True)
    statpop_reduced

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

        statent_filtered[[f'closest_node_{mode}', f'access_egress_{mode}']] = list(zip(*nodes))

    statent_filtered

    # make a sample of statpop cells but inlcude the origin cell in every case
    statent_filtered = pd.concat([
        statent_filtered[statent_filtered['id'] == statent_id],
        statent_filtered[statent_filtered['id'] != statent_id].sample(frac=destinations_sample, random_state=0)
    ])

    tt_matrices = {}

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

    tt_matrices[MODE_TRANSIT]

    # calculate shortest paths to all destinations in networkx
    for mode in [MODE_PRIVATE_CARS, MODE_FOOT, MODE_CYCLING]:

        #print(mode)

        origin_cell = statent_filtered[statent_filtered['id'] == statent_id]
        origin_node = min(origin_cell[f'closest_node_{mode}'])
        access_cost = min(origin_cell[f'access_egress_{mode}'])
        #print(origin_node, access_cost)

        if mode == MODE_PRIVATE_CARS:
            weight = 'matsim_tt_cars'
        else:
            weight = f'cost_{mode}'

        # calculate cost to all other nodes
        dijkstra_path_costs = nx.single_source_dijkstra_path_length(L_modes[mode], origin_node, weight=weight)

        # convert the dijkstra path costs dict to a dataframe like from R5
        tt_matrices[mode] = pd.DataFrame(list(dijkstra_path_costs.items()), columns=['to_node', f'travel_time_{mode}'])
        tt_matrices[mode]['from_node'] = origin_node

    tt_matrices[MODE_PRIVATE_CARS]

    # merge the tt matrix of transit on cell_ids
    statent_filtered_2 = pd.merge(
        statent_filtered, tt_matrices[MODE_TRANSIT],
        left_on='id', right_on='to_cell', how='left'
    )
    statent_filtered_2.drop(columns=['from_cell', 'to_cell'], inplace=True)
    statent_filtered_2[f'travel_time_total_{MODE_TRANSIT}'] = statent_filtered_2[f'travel_time_{MODE_TRANSIT}']
    statent_filtered_2

    # merge the other tt matrices on node ids
    for mode in [MODE_PRIVATE_CARS, MODE_FOOT, MODE_CYCLING]:

        # merge the travel time matrices with the statent dataset
        # such that every destination has the cost of reaching it
        statent_filtered_2 = pd.merge(
            statent_filtered_2, tt_matrices[mode],
            left_on=f'closest_node_{mode}', right_on='to_node', how='left'
        )
        statent_filtered_2.drop(columns=['from_node', 'to_node'], inplace=True)

        if mode == MODE_FOOT:
            speed_factor = walking_speed_kmh / 3.6
        elif mode == MODE_CYCLING:
            speed_factor = cycling_speed_kmh / 3.6
        else:
            speed_factor = 1

        # calculate total travel time by adding access and egress cost and considering speed factors
        statent_filtered_2[f'travel_time_total_{mode}'] = (
                (statent_filtered_2[f'travel_time_{mode}'] / speed_factor)
                + ((statent_filtered_2[f'access_egress_{mode}'] * access_egress_detour_factor) / (
                    walking_speed_kmh / 3.6))
                + ((access_cost * access_egress_detour_factor) / (walking_speed_kmh / 3.6))
        )

    statent_filtered_2

    statent_filtered_2['from_cell'] = statent_id
    statent_filtered_2

    destinations_with_cost = pd.merge(
        statpop_reduced.reset_index(),
        statent_filtered_2.reset_index()[[
            'id', 'VOLLZEITAEQ_TOTAL',
            'travel_time_total_cycling',
            'travel_time_total_foot',
            'travel_time_total_transit',
            'travel_time_total_private_cars',
            'access_egress_cycling',
            'access_egress_private_cars',
            'travel_time_private_cars',
            'closest_node_private_cars',
            'from_cell',
            'geometry'
        ]].rename(columns={'geometry': 'destination_geometry'}),
        how='left', left_on='statent_id', right_on='from_cell')
    destinations_with_cost

    # add base travel time to avoid zero travel time
    for mode in [MODE_CYCLING, MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_FOOT]:
        destinations_with_cost[f'travel_time_total_{mode}'] += base_travel_time

    # calculate behavioral travel cost
    destinations_with_cost[['choice_probabilities', 'behavioral_cost']] = destinations_with_cost.apply(
        lambda row: calculate_behavioral_cost(
            private_cars=row['travel_time_total_private_cars'],
            transit=row['travel_time_total_transit'],
            foot=row['travel_time_total_foot'],
            cycling=row['travel_time_total_cycling'],
            age=row['age']
        ),
        axis=1,
        result_type='expand'
    )

    destinations_with_cost

    # calculate total accessibility contribution using cost function and destination utility
    destinations_with_cost['accessibility_contribution'] = (
            (destinations_with_cost['behavioral_cost'] ** -0.7)
            * destinations_with_cost['VOLLZEITAEQ_TOTAL']
    )

    destinations_with_cost

    # calculate accessibility contribution for every mode separately
    for mode in [MODE_CYCLING, MODE_PRIVATE_CARS, MODE_TRANSIT, MODE_FOOT]:
        destinations_with_cost[f'accessibility_contribution_{mode}'] = (
                (destinations_with_cost[f'travel_time_total_{mode}'] ** -0.7)
                * destinations_with_cost['VOLLZEITAEQ_TOTAL']
        )

    destinations_with_cost

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

    accessibility = gpd.GeoDataFrame(accessibility, crs=2056)
    return accessibility
