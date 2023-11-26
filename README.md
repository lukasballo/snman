# SNMan

**SNMan** (Street Network Manipulator) is a Python toolkit that lets you model urban mobility transitions
by rebuilding street networks according to custom design rules.
For instance, you can use it to model new networks of cycling or transit infrastructure built
by reallocating existing road space from other modes.

It lets you export the updated street networks in various formats, including OSM
that can be used in other tools of your choice for further analysis (e.g., MATSim, SUMO, etc.).

SNMnan is being developed as part of the [E-Bike City](https://ebikecity.baug.ethz.ch/en/) project of
ETH Zurich, Department of Civil and Environmental Engineering.


## Getting Started

  * Install Python 3.9
  * Make sure you have installed the dependencies in the *dependencies* file
  * Pull the *main-v2* branch of this repository
  * Use the Jupyter notebooks under examples, start with *prepare_simplified_street_graph.ipynb*
  * Make sure to change the *data_directory* variable in both files to match the path of your *data_directory*


## Usage

The following Jupyter notebook files are available for a quick start:
  * *prepare_simplified_street_graph.ipynb*: Downloads the network data within
    the specified perimeter from OSM, builds a simplified street graph
    (where each intersection is exacly one node), and enriches the data with other sources, such as
    given public transit routes, elevation, official on-street parking data, traffic volumes, etc.
  * *rebuild_street_graph_multi.ipynb*: Rebuilds the street graph according to a set of rules.

The outputs are generated in *gpkg* format. Use QGIS to view them.


## Changing the process inputs

The following inputs require a geodata starter kit.
Reach out to [Lukas Ballo](https://www.ivt.ethz.ch/personen/profil.lukas-ballo.html to get the newest version.

Input geofiles (under *INPUTS* in *snman_detailed.qgz*).
For an explanation of the columns in each file see
the documentation of each corresponding load function under snman.io, e.g., `snman.io.load_rebuilding_regions()`
  * **rebuilding_regions**: Areas where the streets should be rebuilt
  * **perimeters**: Perimeters for loading OSM data and generating the simplified street network
    (set the active perimeter in *prepare_simplified_street_graph.ipynb* in section *LOAD DATA*)
  * **intersection_polygons**: Manual override of the automatically detected intersection geometries in the simplification
  * **measurement_regions**: Regions where graph measures should be calculated


## Features

SNMan is built on top of [OSMNx](https://osmnx.readthedocs.io/en/stable/), [GeoPandas](https://geopandas.org/), [NetworkX](https://networkx.org/), and [matplotlib](https://matplotlib.org/) and interacts with [OpenStreetMap](https://www.openstreetmap.org/) APIs to:

  * Download street networks
  * Simplify street graphs to obtain exactly one edge per street and one node per intersection (using a process that
    was further developed from OSMNx)
  * Match traffic volumes onto the network
  * Reallocate street space using customizable design rules
  * Automatically arrange a system one-way streets
  * Export as GeoPackage and OSM
  * Calculate metrics of the resulting networks

Planned:

  * Export as [GMNS](https://github.com/zephyr-data-specs/GMNS) format
  * Accessibility analyses and data export for time cartograms


## Contact

For questions or contributions, please contact
[Lukas Ballo](https://www.ivt.ethz.ch/personen/profil.lukas-ballo.html).


## License

SNMan is licensed under the MIT license. OpenStreetMap's open data [license](https://www.openstreetmap.org/copyright/)
requires that derivative works provide proper attribution.


## Citing and further information

  * Ballo, L. and K.W. Axhausen (2023).
    Modeling sustainable mobility futures using an automated process
    of road space reallocation in urban street networks: A case study in Zurich,
    *Transportation Research Board Annual Meeting*, Washington DC, January 2024
