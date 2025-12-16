# SNMan

**SNMan** (Street Network Manipulator) is a Python toolkit that lets you model urban mobility transitions in
a process of rapid experimentation. 
It creates alternative transportation networks through reallocating road space on existing streets. 
For instance, you can use it to model possible networks of cycling or transit infrastructure while ensuring
connectivity and capacity for motorized traffic.

The underlying data is acquired from OpenStreetMap (OSM) and can be optionally enriched by additional data sources,
such as on-street parking, transit routes, or actual street widths.
But SNMan will also work if none of this data is available by extracting as much data as possible from OSM.

For the generation of transport networks, you can define a set of design rules and manual constraints.
By defining smaller rebuilding zones, you can work with street networks of any size, including entire cities.

The updated networks can be exported into standard Geofiles, as well as the OSM format that can be used in
transport modelling tools like MATSim or SUMO.

For details on the underlying methodologies, please refer to the paper at the bottom of this page.

**SNMan** was developed as part of the [E-Bike City](https://ebikecity.baug.ethz.ch/en/) project of
ETH Zurich, Department of Civil and Environmental Engineering.


## Starter data pack and directory structure, download [here](https://polybox.ethz.ch/index.php/s/Uf7LOFrGj3n27oT)

  * **inputs**: Input files such as perimeters, intersection polygons, and
    on-street parking spaces in Zurich
  * **process**: A directory for pre-processed networks
  * **outputs**: Networks after the rebuilding process, based on pre-processed networks
  * **snman_admin_gui.qgz**: A QGIS project file for viewing and editing the input/output files


## Getting Started

For detailed setup instructions, please refer to `documentation/environment_setup.md`. This includes:
  * Installing Python 3.11.13 and setting up the conda environment
  * Installing Cursor IDE and configuring the Python interpreter
  * Installing dependencies (including GDAL via conda)
  * Setting up the `.env` file with your `DATA_DIRECTORY` path
  * Downloading the data starter pack

Once set up, use the Jupyter notebooks under **examples**, see the next section.


## Usage

For a basic setup of your programming environment (Cursor IDE, Git, Miniconda, etc.), please refer to the file
`documentation/environment_setup.md`.

The following Jupyter notebook files are available for a quick start:
  * **prepare_network.ipynb**: Downloads the network data from OSM and prepares it for the rebuilding process
    by running a number of processing steps, including the snman simplification
  * **rebuild_network.ipynb**: Executes the rebuilding process

Please note that the network preparation and the rebuilding process may take several hours to execute,
if you work with networks that cover entire cities.


## Contact

[Lukas Ballo](https://www.ivt.ethz.ch/en/people/profile.lukas-ballo.html)


## License and Citations

SNMan is released under the MIT license. If you use **SNMan** in your work, please cite this paper:

Ballo, L., M. Raubal and K.W. Axhausen (2024)
[Designing an E-Bike City: An automated process for network-wide multimodal road space reallocation](
https://doi.org/10.1016/j.jcmr.2024.100048),
*Journal of Cycling and Micromobility Research*,  **2**, 100048.

