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

For details on the underlying methodologies, please refer to the working paper at the bottom of this page.

**SNMnan** is being developed as part of the [E-Bike City](https://ebikecity.baug.ethz.ch/en/) project of
ETH Zurich, Department of Civil and Environmental Engineering.


## Starter data pack and directory structure, download [here](https://polybox.ethz.ch/index.php/s/Uf7LOFrGj3n27oT)

  * **inputs**: Input files such as perimeters, intersection polygons, and
    on-street parking spaces in Zurich
  * **process**: A directory for pre-processed networks
  * **outputs**: Networks after the rebuilding process, based on pre-processed networks
  * **snman_admin_gui.qgz**: A QGIS project file for viewing and editing the input/output files


## Getting Started


  * Install Python 3.9 and QGIS, we recommend using miniconda
  * Make sure you have installed the dependencies in the **requirements** file.
  * Use the Jupyter notebooks under **examples**, see the next section
  * Install any further missing dependencies


## Usage

For a basic setup of your programming environment (Pycharm, Git, Miniconda, etc.), please refer to the file
`documentation/environment_setup.md`.

The following Jupyter notebook files are available for a quick start:
  * **prepare_network.ipynb**: Downloads the network data from OSM and prepares it for the rebuilding process
    by running a number of processing steps, including the snman simplification
  * **rebuild_network.ipynb**: Executes the rebuilding process

Please note that the network preparation and the rebuilding process may take several hours to execute,
if you work with networks that cover entire cities.


## Contact

**SNMan** is still under development. We put a great effort into documentation. However, if the examples are not
clear enough, please reach out to
[Lukas Ballo](https://www.ivt.ethz.ch/personen/profil.lukas-ballo.html)
for assistance on setting up your working environment.


## License and Citations

SNMan is released under the MIT license. If you use **SNMan** in your work, please cite this working paper:

Ballo, L., M. Raubal and K.W. Axhausen (2024)
[Designing an E-Bike City: An automated process for network-wide multimodal road space reallocation](
https://www.research-collection.ethz.ch/handle/20.500.11850/679713),
*Arbeitsberichte Verkehrs- und Raumplanung*,  **1881**.
