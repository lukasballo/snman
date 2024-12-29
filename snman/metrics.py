import copy

from .constants import *
from . import detours, lane_graph, utils
from . import osmnx_customized as oxc
import pandas as pd
import geopandas as gpd

def metrics_for_measurement_regions_and_synpop(
        measurement_regions,
        synpop,
        include_basic_synpop_metrics=True,
        L=None,
        pois=None,
        synpop_sample=1,
        mode=MODE_PRIVATE_CARS
):

    if synpop_sample < 1:
        synpop['sample_weights'] = synpop['residents'] + synpop['jobs']
        synpop = synpop.sample(
            round(len(synpop) * synpop_sample),
            weights='sample_weights',
            random_state=0
        )

    if not utils.is_in_list(None, [L, pois]):
        synpop['mean_cost'] = synpop['node'].apply(
            lambda node: detours.detour_metrics_for_origin(
                L, node, pois['node'], weight=f'cost_{mode}'
            )['mean_cost']
        )

    agg = {}

    if include_basic_synpop_metrics:
        agg['residents'] = ('residents', 'sum')
        agg['jobs'] = ('stellen_vollzeitaeq', 'sum')

    if not utils.is_in_list(None, [L, pois]):
        agg['mean_of_mean_cost'] = ('mean_cost', 'mean')

    #print('AGG', agg)

    # aggregate synpop to measurement regions
    measurement_regions_synpop = (
        measurement_regions
        .reset_index()
        .sjoin(synpop, predicate='intersects', how='left')
        .groupby('NAME')
        .agg(**agg)
    )

    return measurement_regions_synpop


def metrics_for_measurement_regions_and_lanes(
        L,
        measurement_regions
):

    # Convert the lane graph to GeoDataFrame and aggregate
    lanes = oxc.graph_to_gdfs(L, nodes=False)
    lanes['lane_area'] = lanes.apply(lambda row: row['length'] * row['lane'].width if row['instance'] == 1 else 0, axis=1)
    lanes['lanetype'] = lanes.apply(lambda row: row['lane'].lanetype, axis=1)

    measurement_regions_lanes = (
        measurement_regions
        .reset_index()
        .sjoin(lanes, predicate='intersects', how='left')
        .groupby(['NAME', 'lanetype'])
        .agg(
            lane_area=('lane_area', 'sum')
        )
        .pivot_table(index='NAME', columns='lanetype', values='lane_area')
        .reset_index()
        .set_index('NAME')
    )

    # add a sum of all lane areas
    measurement_regions_lanes['total_lane_area'] = measurement_regions_lanes.sum(axis=1)

    return measurement_regions_lanes


def metrics_for_measurement_regions_and_offstreet_parking(
        measurement_regions,
        offstreet_parking
):

    measurement_regions_parking = (
        measurement_regions
        .reset_index()
        .sjoin(offstreet_parking, predicate='intersects', how='left')
        .groupby(['NAME', 'type'])
        .agg(
            car_parking=('car_parking', 'sum'),
            bicycle_parking=('bicycle_parking', 'sum')
        )
        .pivot_table(index='NAME', columns='type', values=['car_parking', 'bicycle_parking'])
        .reset_index()
        .set_index('NAME')
    )

    measurement_regions_parking.columns = [
        '_'.join(map(str, col)).strip()
        for col
        in measurement_regions_parking.columns.values
    ]

    return measurement_regions_parking


def metrics_for_lane_graph(
        measurement_regions,
        L=None,
        include_measurement_region_metrics=True,
        synpop=None,
        synpop_sample=1,
        pois=None,
        include_basic_synpop_metrics=True,
        include_lane_metrics=True,
        offstreet_parking=None,
        mode=MODE_PRIVATE_CARS,
        buffer=20,
        export_buffered_geometry=False,
        crs=2056
):

    measurement_regions_buffered = copy.deepcopy(measurement_regions)
    measurement_regions_buffered.geometry = measurement_regions.buffer(buffer)
    measurement_regions_buffered['area'] = measurement_regions_buffered.geometry.area

    if include_measurement_region_metrics:
        metrics = measurement_regions_buffered
    else:
        metrics = measurement_regions_buffered[['geometry']]

    # add synpop metrics
    if synpop is not None:
        metrics = pd.merge(
            metrics,
            metrics_for_measurement_regions_and_synpop(
                measurement_regions_buffered,
                synpop,
                synpop_sample=synpop_sample,
                include_basic_synpop_metrics=include_basic_synpop_metrics,
                L=L, pois=pois, mode=mode
            ),
            left_index=True, right_index=True,
            how='left'
        )

    # add lane area metrics
    if L is not None and include_lane_metrics:
        metrics = pd.merge(
            metrics,
            metrics_for_measurement_regions_and_lanes(L, measurement_regions_buffered),
            left_index=True, right_index=True,
            how='left'
        )

    # add offstreet parking metrics
    if offstreet_parking is not None:
        metrics = pd.merge(
            metrics,
            metrics_for_measurement_regions_and_offstreet_parking(measurement_regions_buffered, offstreet_parking),
            left_index=True, right_index=True,
            how='left'
        )

    # replace the buffered geometries with the original ones
    if export_buffered_geometry is False:
        metrics = pd.merge(
            metrics.drop(columns='geometry'),
            measurement_regions[['geometry']],
            left_index=True, right_index=True,
            how='left'
        )

    # Fill all nan values with 0
    metrics.fillna(0, inplace=True)

    return gpd.GeoDataFrame(metrics, geometry='geometry', crs=crs)
