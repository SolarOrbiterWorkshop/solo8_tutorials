# RPW tutorial

This holds the RPW tutorial for the Data Analysis Day at the Solar Orbiter 8th Workshop in Belfast UK.

## Overview

The tutorial is divided in three parts:

1. A short introduction of the RPW instrument science and data products in `about_rpw_data.ipynb`
2. A notebook `get_rpw_data.ipynb` explaining how to retrieve RPW data from SOAR and Sunpy
3. Examples of how RPW data can be used in `rpw_data_tds_tutorial.ipynb` and `rpw_data_tnr_tutorial.ipynb`

## Pre-requisites

Before running the tutorial on your local machine, make sure to install Python package dependencies. 
Here is how to install dependencies using [poetry](https://python-poetry.org/): 

1. Install or upgrade poetry with pip
```
pip install poetry -U
```
2. Then run poetry from the tutorial main directory
```
poetry install
```

NOTES:
* It is also possible to run the tutorial in a [Docker](https://www.docker.com/) container or on mybinder.org --> [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/git/https%3A%2F%2Fgitlab.obspm.fr%2Frpw%2Ftutorials%2Fpython-tnr-tds-tutorial.git/main) (running the tutorial on mybinder.org does not require any action from users concerning installation).
.
* Tutorial relies on Jupyter notebook to work. Visit https://jupyter.org/ for more details.


## Authors

* Antonio Vecchio (antonio.vecchio@obspm.fr) and Xavier Bonnin (xavier.bonnin@obspm.fr) with the support of the SolO/RPW team.

## License

RPW tutorial is released under MIT license.

