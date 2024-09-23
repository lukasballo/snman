This file will guide you through the process of setting up your working environment before working with **SNMan**.
The instructions apply for Windows users.

We will use the following setup:
- Pycharm IDE
- Miniconda3
- Git
- Java
- The data starter pack


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
Remember the path of you mininiconda3 installation, something like `C:\Users\ukas\miniconda3\`

Open Command Prompt, type 'conda' and press enter. If you get a message like 'not recognized', then you need to
add miniconda to the system path:

1. Open Start and type 'path'
2. Click on 'Edit system environment variables'
3. Click on 'Environment variables'
4. Under system variables (bottom half), click on 'path' and click 'edit...'
5. Click on 'New' and add the conda path, usually something like `C:\Users\lukas\miniconda3\Scripts`
6. Reboot

Now, type the following. This will create a new environment in miniconda and activate it:

```
conda create -n snman python=3.9
activate snman
```


# Pycharm #

Download and install Pycharm: https://www.jetbrains.com/pycharm/ The free community version is sufficient but
your institution may also have access to the better, professional version.

Open Pycharm and create a new project in the repository from above (File > New Project). Select 'Pure Python'.

Go to File > Settings > Project:snman > Python interpreter > snman > Add interpreter >
Add local interpreter > Conda Environment. Then, for conda executable, set the path of your miniconda3 installation,
typically `C:\Users\lukas\miniconda3\Scripts\conda.exe`

Click on 'Load Environments', select 'Existing environment' and 'snman'.

Now, open the terminal in Pycharm (one of the buttons bottom left) and install the dependencies in the
**requirements** file. For each dependency,
first try to install it using conda, such as the following example:

```
conda install geopandas
```

Some packages are not available on conda (e.g., osmnx, r5py). For these use pip:

```
pip install osmnx
```

Now, try to run the Jupyter Notebook `examples/prepare_network.ipynb`

If you get an error about missing packages, install them using the steps above.


# Java #

If you get an error message like
`JVMNotFoundException: No JVM shared library found...`, download and install Java and set the `JAVA_HOME` path using the following steps

1. Search for your java installation path, it will be typically something like
   `C:\Program Files (x86)\Common Files\Java\jre1.8.0_421\bin\`.
2. Open Start and type 'path'
3. Click on 'Edit system environment variables'
4. Click on 'Environment variables'
5. Under system variables (bottom half), click on 'add...'
6. name it `JAVA_HOME` and paste the path
7. Reboot