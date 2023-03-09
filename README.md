# SNMan

**SNMan** is an under-development Python toolkit that lets you analyze street networks and model radical urban mobility transitions. It lets you export the updated street networks in OSM format that can be used in a transport modelling tool of your choice (e.g., MATSim, SUMO, etc.).

A detailed description is underway in an upcoming research paper:

  * Ballo, L. 2023. "SNMan: Modeling radical urban mobility transitions"
    *23rd Swiss Transport Research Conference*, Ascona, May 2023

SNMnan is being developed as part of the [E-Bike City](https://ebikecity.baug.ethz.ch/en/) project of
ETH Zurich, Department of Civil and Environmental Engineering.


## Getting Started

  * Install Python 3.9
  * Make sure you have installed the following dependencies: geopandas, osmnx, shapely, statistics, itertools, gdal, rasterio
  * Pull the *main* branch of this repository
  * [Download](https://polybox.ethz.ch/index.php/s/2yjdcNX1kJmgw8W) the *data_directory* starter kit
    and save it on your local machine
  * Start with *prepare_simplified_street_graph.ipynb* and *rebuild_street_graph.ipynb* for a complete process
    of network simplification, enrichment, and rebuilding into a system of one-way streets. 
    Change the *data_directory* variable in both files to match the path of your *data_directory*
  * A simplified and enriched network of Zurich is provided in the starter kit. 
    Open *cached_network.qgz* in [QGIS](https://qgis.org/) to view it on a map.
  * Run *rebuild_street_graph.ipynb*
    to generate a network with rebuilt road space.
  * Open *ebc_preview.qgz* in QGIS to view the resulting network


## Changing the process inputs

Input geofiles (under *INPUTS* in *cached_network.qgz* and *ebc_preview.qgz*).
For an explanation of the columns in each file see
the documentation of each corresponding load function under snman.io, e.g., `snman.io.load_rebuilding_regions()`
  * **rebuilding_regions**: Areas where the streets should be rebuilt
  * **perimeters**: Perimeters for loading OSM data and generating the simplified street network
    (set the active perimeter in *prepare_simplified_street_graph.ipynb* in section *LOAD DATA*)
  * **intersection_polygons**: Manual override of the automatically detected intersection geometries in the simplification
  * **regions**: Local override of the network simplification settings


## Features

OSMnx is built on top of [OSMNx](https://osmnx.readthedocs.io/en/stable/), [GeoPandas](https://geopandas.org/), [NetworkX](https://networkx.org/), and [matplotlib](https://matplotlib.org/) and interacts with [OpenStreetMap](https://www.openstreetmap.org/) APIs to:

  * Download street networks
  * Simplify street graphs to obtain exactly one edge per street and one node per intersection (using a process that
    was further developed from OSMNx)
  * Match traffic volumes onto the network
  * Reallocate street space using customizable heuristics
  * Automatically arrange a system one-way streets
  * Export as GeoPackage and OSM

Planned:

  * Export as [GMNS](https://github.com/zephyr-data-specs/GMNS) format
  * Accessibility analyses and data export for time cartograms


## Contact

For questions or contributions, please contact
[Lukas Ballo](https://www.ivt.ethz.ch/personen/profil.lukas-ballo.html).


## License

SNMan is licensed under the MIT license. OpenStreetMap's open data [license](https://www.openstreetmap.org/copyright/)
requires that derivative works provide proper attribution.

