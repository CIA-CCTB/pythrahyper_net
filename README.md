# pythrahyper_net
Code and data for the paper *"Biological network growth in complex environments - a computational framework"* by T. Paul and P. Kollmannsberger (2020) - https://biorxiv.org/cgi/content/short/2020.06.01.127407v1

Please have a look at the notebook [Introduction.ipynb](https://github.com/CIA-CCTB/pythrahyper_net/blob/master/Introduction.ipynb), or try it directly here:  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/CIA-CCTB/pythrahyper_net/blob/master/Colab/Introduction_Colab.ipynb) (requires Google account)

The following Jupyter notebooks reproduce the simulations shown in Figure 6 in the paper:

- [Multicellular-Network.ipynb](https://github.com/CIA-CCTB/pythrahyper_net/blob/master/Multicellular-Network.ipynb) - simulation of network growth between layers of cells
- [Multi-Simulation-Setup.ipynb](https://github.com/CIA-CCTB/pythrahyper_net/blob/master/Multi-Simulation-Setup.ipynb) - generate configuration and batch files for parameter scan
- [Multi-Simulation-Analysis.ipynb](https://github.com/CIA-CCTB/pythrahyper_net/blob/master/Multi-Simulation-Analysis.ipynb) - analyze results and generate plots for parameter scan

# Instructions

The framework is written in python using numpy and the multiprocessing module, and has been tested under Linux and MacOS. To run the example notebooks, first download or clone this repository, and then follow the instructions below.

## 1) Using conda

The easiest way to install the required python packages is by using conda. Creating a new environment with this command will install all dependencies:

`conda create --name pythra python=3.7 pyqt=5 numpy scipy tifffile jupyter networkx matplotlib`

Then change into the new environment using `conda activate pythra`, and start a Jupyter notebook server in the `pythrahyper_net` directory to access the notebooks.

### Mayavi visualization in the browser:

To get interactive mayavi visualizations inside the browser, first install mayavi and ipyevents:

`conda install -c anaconda mayavi`

`conda install -c conda-forge ipyevents`

Next, install and activate the required extension:

`jupyter nbextension install --py mayavi --user`

`jupyter nbextension enable --py mayavi --user`

If you get missing symbol errors upon importing `mlab`, try this:

`conda install -c conda-force "libnetcdf=4.6.2"`

### Interactive Matplotlib plots:

The matplotlib plots can be made interactive using these modules:

`conda install -c conda-forge ipympl widgetsnbextension`

## 2) Using Singularity container

The second possibility is to run the framework inside a Singularity container. A container image can be created using the included definition file:

`sudo singularity build pythra.simg Singularity.def`

After successful build, you can e.g. start a Jupyter notebook server inside the container:

`singularity exec pythra.simg jupyter notebook`

Then copy and paste the server URL into a web browser running outside of the container to access the notebooks.


