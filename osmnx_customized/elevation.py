"""Get node elevations and calculate edge grades."""

import multiprocessing as mp
import time
from hashlib import sha1
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import geopandas as gpd
import requests

from . import downloader
from . import utils
from . import utils_graph

# rasterio and gdal are optional dependencies for raster querying
try:
    import rasterio
    from osgeo import gdal
except ImportError:  # pragma: no cover
    rasterio = gdal = None


def _query_raster(nodes, filepath, band):
    """
    Query a raster for values at coordinates in a DataFrame's x/y columns.

    Parameters
    ----------
    nodes : pandas.DataFrame
        DataFrame indexed by node ID and with two columns: x and y
    filepath : string or pathlib.Path
        path to the raster file or VRT to query
    band : int
        which raster band to query

    Returns
    -------
    nodes_values : zip
        zipped node IDs and corresponding raster values
    """
    # must open raster file here: cannot pickle it to pass in multiprocessing
    with rasterio.open(filepath) as raster:
        values = np.array(tuple(raster.sample(nodes.values, band)), dtype=float).squeeze()
        values[values == raster.nodata] = np.nan
        return zip(nodes.index, values)


def add_node_elevations_raster(G, filepath, band=1, cpus=None, graph_crs=2056, raster_crs=4326):
    """
    Add `elevation` attribute to each node from local raster file(s).

    If `filepath` is a list of paths, this will generate a virtual raster
    composed of the files at those paths as an intermediate step.

    See also the `add_edge_grades` function.

    snman.oxc: added custom raster_crs

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph, in same CRS as raster
    filepath : string or pathlib.Path or list of strings/Paths
        path (or list of paths) to the raster file(s) to query
    band : int
        which raster band to query
    cpus : int
        how many CPU cores to use; if None, use all available
    graph_crs : int
        CRS of the graph geometries
    raster_crs : int
        CRS of the raster file

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with node elevation attributes
    """
    if rasterio is None or gdal is None:  # pragma: no cover
        raise ImportError("gdal and rasterio must be installed to query raster files")

    if cpus is None:
        cpus = mp.cpu_count()
    cpus = min(cpus, mp.cpu_count())
    utils.log(f"Attaching elevations with {cpus} CPUs...")

    # if a list of filepaths is passed, compose them all as a virtual raster
    # use the sha1 hash of the filepaths list as the vrt filename
    if not isinstance(filepath, (str, Path)):
        filepaths = [str(p) for p in filepath]
        sha = sha1(str(filepaths).encode("utf-8")).hexdigest()
        filepath = f"./.osmnx_{sha}.vrt"
        gdal.BuildVRT(filepath, filepaths).FlushCache()

    # convert the x and y coordinates of the nodes to match the raster image and then extract them
    nodes = utils_graph.graph_to_gdfs(G, edges=False, node_geometry=False)
    nodes = gpd.GeoDataFrame(
        nodes, geometry=gpd.points_from_xy(nodes.x, nodes.y), crs=graph_crs
    )
    nodes = nodes.to_crs(raster_crs)
    nodes.x = nodes.geometry.x
    nodes.y = nodes.geometry.y
    nodes = nodes[["x", "y"]]
    print(nodes)

    if cpus == 1:
        elevs = dict(_query_raster(nodes, filepath, band))
    else:
        # divide nodes into equal-sized chunks for multiprocessing
        size = int(np.ceil(len(nodes) / cpus))
        args = ((nodes.iloc[i : i + size], filepath, band) for i in range(0, len(nodes), size))
        pool = mp.Pool(cpus)
        sma = pool.starmap_async(_query_raster, args)
        results = sma.get()
        pool.close()
        pool.join()
        elevs = {k: v for kv in results for k, v in kv}

    assert len(G) == len(elevs)
    nx.set_node_attributes(G, elevs, name="elevation")
    utils.log("Added elevation data from raster to all nodes.")
    return G


def add_node_elevations_google(
    G, api_key, max_locations_per_batch=350, pause_duration=0, precision=3
):  # pragma: no cover
    """
    Add `elevation` (meters) attribute to each node using a web service.

    This uses the Google Maps Elevation API and requires an API key. For a
    free, local alternative, see the `add_node_elevations_raster` function.
    See also the `add_edge_grades` function.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph
    api_key : string
        a Google Maps Elevation API key
    max_locations_per_batch : int
        max number of coordinate pairs to submit in each API call (if this is
        too high, the server will reject the request because its character
        limit exceeds the max allowed)
    pause_duration : float
        time to pause between API calls, which can be increased if you get
        rate limited
    precision : int
        decimal precision to round elevation values

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with node elevation attributes
    """
    # elevation API endpoint ready for use
    url_template = "https://maps.googleapis.com/maps/api/elevation/json?locations={}&key={}"

    # make a pandas series of all the nodes' coordinates as 'lat,lng'
    # round coordinates to 5 decimal places (approx 1 meter) to be able to fit
    # in more locations per API call
    node_points = pd.Series(
        {node: f'{data["y"]:.5f},{data["x"]:.5f}' for node, data in G.nodes(data=True)}
    )
    n_calls = int(np.ceil(len(node_points) / max_locations_per_batch))
    utils.log(f"Requesting node elevations from the API in {n_calls} calls")

    # break the series of coordinates into chunks of size max_locations_per_batch
    # API format is locations=lat,lng|lat,lng|lat,lng|lat,lng...
    results = []
    for i in range(0, len(node_points), max_locations_per_batch):
        chunk = node_points.iloc[i : i + max_locations_per_batch]
        locations = "|".join(chunk)
        url = url_template.format(locations, api_key)

        # check if this request is already in the cache (if global use_cache=True)
        cached_response_json = downloader._retrieve_from_cache(url)
        if cached_response_json is not None:
            response_json = cached_response_json
        else:
            try:
                # request the elevations from the API
                utils.log(f"Requesting node elevations: {url}")
                time.sleep(pause_duration)
                response = requests.get(url)
                response_json = response.json()
                downloader._save_to_cache(url, response_json, response.status_code)
            except Exception as e:
                utils.log(e)
                utils.log(f"Server responded with {response.status_code}: {response.reason}")

        # append these elevation results to the list of all results
        results.extend(response_json["results"])

    # sanity check that all our vectors have the same number of elements
    if not (len(results) == len(G) == len(node_points)):
        raise Exception(
            f"Graph has {len(G)} nodes but we received {len(results)} results from elevation API"
        )
    else:
        utils.log(
            f"Graph has {len(G)} nodes and we received {len(results)} results from elevation API"
        )

    # add elevation as an attribute to the nodes
    df = pd.DataFrame(node_points, columns=["node_points"])
    df["elevation"] = [result["elevation"] for result in results]
    df["elevation"] = df["elevation"].round(precision)
    nx.set_node_attributes(G, name="elevation", values=df["elevation"].to_dict())
    utils.log("Added elevation data from Google to all nodes.")

    return G


def add_edge_grades(G, add_absolute=True, precision=3):
    """
    Add `grade` attribute to each graph edge.

    Vectorized function to calculate the directed grade (ie, rise over run)
    for each edge in the graph and add it to the edge as an attribute. Nodes
    must already have `elevation` attributes to use this function.

    See also the `add_node_elevations_raster` and `add_node_elevations_google`
    functions.

    Parameters
    ----------
    G : networkx.MultiDiGraph
        input graph with `elevation` node attribute
    add_absolute : bool
        if True, also add absolute value of grade as `grade_abs` attribute
    precision : int
        decimal precision to round grade values

    Returns
    -------
    G : networkx.MultiDiGraph
        graph with edge `grade` (and optionally `grade_abs`) attributes
    """
    elev_lookup = G.nodes(data="elevation")
    u, v, k, lengths = zip(*G.edges(keys=True, data="length"))
    uvk = tuple(zip(u, v, k))

    # calculate edges' elevation changes from u to v then divide by lengths
    elevs = np.array([(elev_lookup[u], elev_lookup[v]) for u, v, k in uvk])
    grades = ((elevs[:, 1] - elevs[:, 0]) / np.array(lengths)).round(precision)
    nx.set_edge_attributes(G, dict(zip(uvk, grades)), name="grade")

    # optionally add grade absolute value to the edge attributes
    if add_absolute:
        nx.set_edge_attributes(G, dict(zip(uvk, np.abs(grades))), name="grade_abs")

    utils.log("Added grade attributes to all edges.")
    return G
