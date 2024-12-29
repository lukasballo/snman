import numpy as np
import geopandas as gpd
import shapely as shp
from .constants import *


def fit_edge(
        G, u, v, k, street_polygons_gdf,
        default_width_key=KEY_LANES_DESCRIPTION + '_width_total_m',
        length_to_remove_at_ends=10,
        replace_geometry=False
):
    """
    Fits the edge into the street polygons

    Parameters
    ----------
    G: nx.MultiDiGraph
        street graph
    u: int
    v: int
    k: int
    street_polygons_gdf: gpd.GeoDataFrame
        a geodataframe with street polygons
    default_width_key: str
        which key should be used for width if the fitting is unsuccessful
    replace_geometry: bool
        if true, the edge geometries will be replaced by the matched street axis.
        for debugging purposes only.

    Returns
    -------
    """

    G.edges[(u, v, k)]['width'] = G.edges[(u, v, k)][default_width_key]
    line = G.edges[(u, v, k)]['geometry']

    if G.edges[(u, v, k)].get('layer', 0) != 0:
        G.edges[(u, v, k)]['_fitting_status'] = 'not_on_ground'
        return

    try:

        multipolygon = street_polygons_gdf[street_polygons_gdf.intersects(line.buffer(50))].unary_union

        offset_line_1 = line.parallel_offset(20)
        offset_line_2 = line.parallel_offset(-20)

        line_length_rounded = round(line.length)
        step_length = 10
        interpolation_steps = np.array(
            range(0, line_length_rounded, step_length)
        )

        # remove steps at the beginning and at the end
        interpolation_steps = interpolation_steps[length_to_remove_at_ends < interpolation_steps]
        interpolation_steps = interpolation_steps[interpolation_steps < line_length_rounded - length_to_remove_at_ends]

        cross_sections = [
            shp.LineString([
                offset_line_1.interpolate(position, normalized=True),
                offset_line_2.interpolate(position, normalized=True)
            ])
            for position
            in interpolation_steps / line_length_rounded
        ]

        def intersect_cross_section(cross_section):
            a = cross_section.intersection(multipolygon)
            return a

        cross_sections = [intersect_cross_section(cross_section) for cross_section in cross_sections]
        #print(cross_sections)
        cross_sections = list(
            filter(lambda x: type(x) in (shp.LineString, shp.MultiLineString) and not x.is_empty, cross_sections))

        #if len(cross_sections) < 4:
        #    return

        cross_sections = cross_sections[1:-1]

        def total_length(multiline):
            if type(multiline) == shp.MultiLineString:
                return sum([a.length for a in multiline.geoms])
            elif type(multiline == shp.LineString):
                return multiline.length
            else:
                return 0

        widths = [total_length(cross_section) for cross_section in cross_sections]
        axis_points = [cross_section.centroid for cross_section in cross_sections]

        # left_side_points = [cross_section.interpolate(0, normalized=True) for cross_section in cross_sections]
        # right_side_points = [cross_section.interpolate(1, normalized=True) for cross_section in cross_sections]
        # fitted_polygon = shp.Polygon(left_side_points + list(reversed(right_side_points)))

        if replace_geometry:
            axis = shp.LineString(axis_points)
            G.edges[(u, v, k)]['geometry'] = axis

        G.edges[(u, v, k)]['width'] = np.percentile(widths, 50)
        G.edges[(u, v, k)]['widths_along'] = str(widths)
        G.edges[(u, v, k)]['_fitting_status'] = 'successful'

    except:
        #print('Fitting error', u, v, k)
        pass


def fit_edges(G, street_polygons_gdf, verbose=False, **kwargs):
    """
    Fit the edges into street polygons and replace the width accordingly.

    Parameters
    ----------
    G: nx.MultiDiGraph
    street_polygons_gdf: gpd.GeoDataFrame

    Returns
    -------
    """

    for uvk, data in G.edges.items():
        if verbose:
            print(uvk)
        fit_edge(G, *uvk, street_polygons_gdf, **kwargs)