This file will guide you through the process of setting up your working environment before working with **SNMan**.
The instructions are primarily written for Windows users, but can be adapted for macOS and Linux (e.g., use `conda activate snman` instead of `activate snman` on Unix systems).

We will use the following setup:
- Cursor IDE
- Miniconda3
- Git
- Java (optional, required only if using r5py for transit accessibility calculations)
- QGIS (optional, for viewing and editing geospatial data)
- The data starter pack (download [here](https://polybox.ethz.ch/index.php/s/Uf7LOFrGj3n27oT))


# Git #

Download and install git: https://git-scm.com/downloads

Create a new empty folder for your snman code, such as `C:\Users\lukas\git\snman`

Open command prompt and type the following. It will create a local git repository and pull the code from GitHub. 

```
cd C:\Users\lukas\git\snman
git init
git remote add origin https://github.com/lukasballo/snman
git pull origin main
```


# Miniconda3 #

Download and install Miniconda3: https://docs.anaconda.com/miniconda/#miniconda-latest-installer-links
Remember the path of your miniconda3 installation, something like `C:\Users\lukas\miniconda3\`

Open Command Prompt, type 'conda' and press enter. If you get a message like 'not recognized', then you need to
add miniconda to the system path:

1. Open Start and type 'path'
2. Click on 'Edit system environment variables'
3. Click on 'Environment variables'
4. Under system variables (bottom half), click on 'path' and click 'edit...'
5. Click on 'New' and add the conda path, usually something like `C:\Users\lukas\miniconda3\Scripts`
6. Reboot

Now, type the following. This will create a new environment in miniconda with Python 3.11.13 and activate it:

**Windows:**
```
conda create -n snman python=3.11.13
activate snman
```

**macOS/Linux:**
```
conda create -n snman python=3.11.13
conda activate snman
```

**Note:** SNMan requires Python 3.11.13. Make sure to use this exact version for compatibility.


# Cursor IDE #

Download and install Cursor: https://cursor.sh/

Open Cursor and open the repository folder from above (File > Open Folder).

To configure the Python interpreter:
1. Open the Command Palette (Ctrl+Shift+P or Cmd+Shift+P)
2. Type "Python: Select Interpreter"
3. Select the conda environment, typically located at `C:\Users\lukas\miniconda3\envs\snman\python.exe`


# Install Dependencies #

Now, open the terminal in Cursor and install the required packages. First, install gdal using conda:

```
conda install -c conda-forge gdal
```

Create a `.env` file in the project root based on `.env.example` and set your `DATA_DIRECTORY` path.

Then, automatically install all dependencies with exact versions from the **requirements** file:

```
pip install -r requirements
```

The **requirements** file contains all required Python packages with their exact pinned versions to ensure reproducibility. The Python version (3.11.13) is also specified at the top of the requirements file.

**Note:** GDAL should be installed via conda (as shown above) rather than pip, as it has complex system dependencies that are better handled by conda.

Now, try to run the Jupyter Notebook `examples/prepare_network.ipynb`.
If you encounter any missing packages (which should not happen if the requirements file is complete), you can install them manually:

```
pip install <package name>
```


# Java (Optional) #

Java is only required if you plan to use `r5py` for transit accessibility calculations. The `r5py` import is optional in SNMan, so Java is not strictly necessary for basic functionality.

If you need Java and get an error message like `JVMNotFoundException: No JVM shared library found...` or `CalledProcessError` related to Java, download and install Java (JDK 11 or later) and set the `JAVA_HOME` environment variable.

**Windows:**
1. Download and install Java JDK from https://adoptium.net/ or your preferred Java distribution
2. Find your Java installation path, typically something like `C:\Program Files\Microsoft\jdk-XX.X.X-hotspot\` or `C:\Program Files\Java\jdk-XX`
3. Open Start and type 'path'
4. Click on 'Edit system environment variables'
5. Click on 'Environment variables'
6. Under system variables (bottom half), click on 'add...'
7. Name it `JAVA_HOME` and paste the path to the JDK directory (without the `bin` subdirectory)
8. Reboot

**macOS/Linux:**
1. Install Java using your system package manager or download from https://adoptium.net/
2. Set `JAVA_HOME` in your shell configuration file (e.g., `~/.zshrc` or `~/.bashrc`):
   ```
   export JAVA_HOME=/path/to/java
   ```
3. Reload your shell configuration or restart your terminal


# Data Starter Pack #

Download the starter data pack [here](https://polybox.ethz.ch/index.php/s/Uf7LOFrGj3n27oT). This pack contains:

  * **inputs**: Input files such as perimeters, intersection polygons, and
    on-street parking spaces in Zurich
  * **process**: A directory for pre-processed networks
  * **outputs**: Networks after the rebuilding process, based on pre-processed networks
  * **snman_admin_gui.qgz**: A QGIS project file for viewing and editing the input/output files

After downloading, extract the data pack and set the `DATA_DIRECTORY` path in your `.env` file to point to the extracted directory.