from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
import geopandas as gpd
from .constants import *
import datetime
import r5py


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
        transport_network_default,
        transport_network_car,
        transport_network_cycling,
        statpop_with_statent_ids,
        statent,
        statent_id,
        distance_limit=10000,
        base_travel_time=0.1,
        destinations_sample=1,
        population_sample=1
):
    """
    Calculates accessibility for every resident associated with a given cell in the statent dataset.

    Parameters
    ----------
    transport_network_default: r5py.TransportNetwork
        transport network without travel time adjustments
    transport_network_car: r5py.TransportNetwork
        transport network for cars, with maxspeeds adjusted to reflect the congested state from MATSim
    transport_network_cycling: r5py.TransportNetwork
        transport network for cyclists, with all usable links encoded as highway=primary and maxspeeds such that
        the resulting link travel times reflect the VoD-adjusted lengths
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
        travrel time in minutes to be added to all travel timed, a small number is useful to prevent 0 min. travel times

    Returns
    -------
    gpd.GeoDataFrame
    """

    # filter statpop for residents within the given statent cell
    statpop_filtered = statpop_with_statent_ids.query(f'statent_id == {statent_id}')
    statpop_filtered = statpop_filtered.sample(frac=population_sample, random_state=0)
    if len(statpop_filtered) == 0:
        return None

    # filter statent for cells within a distance limit
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

    statent_filtered = pd.concat([
        statent_filtered[statent_filtered['id'] == statent_id],
        statent_filtered[statent_filtered['id'] != statent_id].sample(frac=destinations_sample, random_state=0)
        ])

    # calculate travel time matrices
    tt_computers = {}
    common_settings = {
        'origins': origins,
        'destinations': statent_filtered,
        'departure': datetime.datetime(2024, 2, 22, 8, 30),
        'snap_to_network': True
    }

    tt_computers[MODE_FOOT] = r5py.TravelTimeMatrixComputer(
        transport_network_default,
        transport_modes=[r5py.TransportMode.WALK],
        **common_settings
    )

    # cycling trips are routed using the "car" mode, on a network with following specialties:
    # - all links usable by cyclists are encoded as highway=primary
    # - maxspeeds are calculated such that the resulting travel times reflect the VoD-adjusted lengths
    tt_computers[MODE_CYCLING] = r5py.TravelTimeMatrixComputer(
        transport_network_cycling,
        transport_modes=[r5py.TransportMode.CAR],
        **common_settings
    )

    # the "car" network has the following adjustments:
    # - maxspeeds calculated such that the resulting travel times correspond to the congested state from MATSim
    tt_computers[MODE_PRIVATE_CARS] = r5py.TravelTimeMatrixComputer(
        transport_network_car,
        transport_modes=[r5py.TransportMode.CAR],
        **common_settings
    )

    tt_computers[MODE_TRANSIT] = r5py.TravelTimeMatrixComputer(
        transport_network_default,
        transport_modes=[r5py.TransportMode.TRANSIT],
        **common_settings
    )

    def tt_calculation(tt_computer, mode):
        tt_matrix = tt_computer.compute_travel_times()
        # tt_computer.compute_travel_details()
        tt_matrix['mode'] = mode
        return tt_matrix

    # calculate a travel time matrix for each mode
    tt_matrices = {mode: tt_calculation(tt_computer, mode) for mode, tt_computer in tt_computers.items()}

    # merge all travel time matrices into one
    tt_matrix = pd.concat(tt_matrices.values())

    # add base travel time to avoid 0 values
    tt_matrix['travel_time'] = tt_matrix['travel_time'] + base_travel_time

    # make a light version of statpop
    statpop_reduced = statpop_filtered.reset_index()[['record', 'age', 'statent_id', 'geometry']]
    statpop_reduced.set_index('record', inplace=True)

    accessibility_costs = pd.merge(statpop_reduced.reset_index(), tt_matrix, how='left', left_on='statent_id',
                                   right_on='from_id')
    accessibility_costs = accessibility_costs[accessibility_costs['travel_time'].notna()]

    # group accessibility by person, providing a list of travel times to each statent destination

    ag = accessibility_costs.reset_index().groupby(['record', 'to_id']).agg({
        'travel_time': list,
        'mode': list,
        # 'sex': 'first',
        'age': 'first',
        # 'maritalstatus': 'first',
        # 'residencepermit': 'first',
        # 'residentpermit': 'first'
    })

    # add destination information from statent
    ag = pd.merge(ag.reset_index(), statent.reset_index()[['id', 'VOLLZEITAEQ_TOTAL']], left_on='to_id', right_on='id')

    # create a dictionary with travel options for each person and statent destination
    ag['travel_options'] = ag.apply(lambda row: dict(zip(row['mode'], row['travel_time'])), axis=1)

    # run the mode choice model to calculate each person's generalized cost to each statent destination
    ag[['choice_probabilities', 'behavioral_cost']] = ag.apply(
        lambda row: calculate_behavioral_cost(**row['travel_options'], age=row['age']),
        axis=1,
        result_type='expand'
    )

    # apply the accessibility cost function
    ag['accessibility_contribution'] = ag['behavioral_cost'] ** -0.7

    # for each person, sum the accessibility contributions across all destinations
    accessibility = ag.groupby('record').agg({
        'accessibility_contribution': 'sum'
    })

    accessibility.rename(columns={'accessibility_contribution': 'accessibility'}, inplace=True)
    accessibility['accessibility'] = accessibility['accessibility'] / destinations_sample

    accessibility = accessibility.join(statpop_reduced)
    return gpd.GeoDataFrame(accessibility, crs=2056)


from joblib import Parallel, delayed
import time

# Define a sample function to be executed
def my_function(args, x):
    return calculate_accessibility_for_statent_cell(*args, x)

def parallel(args, inputs):
    results = Parallel(n_jobs=3)(delayed(my_function)(i) for i in inputs)

