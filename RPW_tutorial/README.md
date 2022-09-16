# RPW tutorial

This holds the RPW tutorial for the Data Analysis Day at the Solar Orbiter 8th Workshop in Belfast UK.

## Overview

The tutorial is divided in three parts:

1. A short introduction of the RPW instrument science and data products in `about_rpw_data.ipynb`
2. A notebook `get_rpw_data.ipynb` explaining how to retrieve RPW data from SOAR and Sunpy
3. Examples of how RPW data can be used in `rpw_data_tds_tutorial.ipynb` and `rpw_data_tnr_tutorial.ipynb`

## User guide

### Running the tutorial on a local machine

Before running the tutorial on your local machine, make sure to have the following software installed:
- Python 3 and pip (tested on Python 3.8)
- gfortran
- wget
- tar

To install [NASA CDF software](https://cdf.gsfc.nasa.gov/), enter:

```/bin/bash install_cdf.sh```

To install Python package dependencies, enter:

```/bin/bash install.sh <venv>```

Where `<venv>` is the name of the [virtual environment](https://docs.python.org/3/tutorial/venv.html) to use to run the tutorial (default is `<venv>=rpw_tuto`)
The virtualenv source files will be stored in `~/.virtualenvs/<venv>`.

Once installed, you can start jupyter using command:

```/bin/bash start_tuto.sh <venv>```

When running a notebook, make sure to use the right kernel (`<venv>`). 
Kernel can be changed from the top menu. Visit https://jupyter.org/ for more details.

NOTES:
* It is also possible to build and run the notebook in a Docker container. Visit [Docker](https://www.docker.com/) to learn how to use Docker.

### Running the tutorial on Binder

You can try to run the tutorial on Blinder, just click on [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fgitlab.obspm.fr%2Frpw%2Ftutorials%2Fpython-tnr-tds-tutorial.git/main)

NOTE: the installation is automatically performed by Blinder using Dockerfile. 

## Authors

* Antonio Vecchio (antonio.vecchio@obspm.fr) and Xavier Bonnin (xavier.bonnin@obspm.fr) with the support of the SolO/RPW team.

## License

RPW tutorial is released under MIT license.

