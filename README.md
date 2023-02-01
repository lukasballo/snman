# SNMan

**SNMan** is an under-development Python toolkit that lets you analyze street networks and model radical urban mobility transitions. It lets you export the updated street networks in OSM format that can be used in a transport modelling tool of your choice (e.g., MATSim, SUMO, etc.).

A detailed description is underway in an upcoming research paper:

  * Ballo, L. 2023. "SNMan: Modeling radical urban mobility transitions." *23rd Swiss Transport Research Conference*, Ascona, May 2023

SNMnan is being developed as part of the [E-Bike City](https://ebikecity.baug.ethz.ch/en/) project of ETH Zurich, Department of Civil and Environmental Engineering:


## Getting Started

  * Install Python 3.9
  * Pull the *main* branch of this repository
  * Obtain the GeoData starter kit from [Lukas Ballo](https://www.ivt.ethz.ch/personen/profil.lukas-ballo.html) (it will be published here later)
  * Run *main_pipeline.py*

## Features

OSMnx is built on top of [OSMNx](https://osmnx.readthedocs.io/en/stable/), [GeoPandas](https://geopandas.org/), [NetworkX](https://networkx.org/), and [matplotlib](https://matplotlib.org/) and interacts with [OpenStreetMap](https://www.openstreetmap.org/) APIs to:

  * Download street networks
  * Simplify street graphs to obtain exactly one edge per street and one node per intersection (using a process that was further developed from OSMNx)
  * Reallocate street space using customizable heuristics
  * Automatically arrange a system one-way streets
  * Export as GeoPackage and OSM

Planned:

  * Matching traffic volumes on the street network
  * Export as [GMNS](https://github.com/zephyr-data-specs/GMNS) format

## Roadmap

For the feature roadmap, see here: https://trello.com/b/XfO5BLRg/snman

## Contact

For questions or contributions, please contact
[Lukas Ballo](https://www.ivt.ethz.ch/personen/profil.lukas-ballo.html).

## License

SNMan is licensed under the MIT license. OpenStreetMap's open data [license](https://www.openstreetmap.org/copyright/) requires that derivative works provide proper attribution.

