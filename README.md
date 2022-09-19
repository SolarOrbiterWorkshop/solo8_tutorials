# Solar Orbiter 8 Workshop [Data Analysis Tutorial Day](https://blogs.qub.ac.uk/so8belfast/data-analysis-workshop-16-september-2022-2/) 16 Sept 2022


<div>
<img src="./images/Solar_Orbiter_reaches_first_perihelion_pillars.jpeg" width="300" align="left"/>
</div>

This is a repository to hold all the tutorial notebooks for the data analysis day of the Solar Orbiter 8 Workshop. 

The notebooks for each instrument is in their repective folders.

*NOTE*:  For the **EUI**, **Metis**,  **PHI** and  **Solo-HI** tutorials accompanying data will need to be downloaded to run the tutorial notebooks. The links to download these data are within the README files in the instrument directories.


## JHelioviewer tutorial

JHelioviewer can be downloaded from here: https://www.jhelioviewer.org/

The movies and he JHelioviewer states that were shown in this tutorial can be accessed here:
 * ftp://omaftp.oma.be/dist/astro/Verbeeck.F/JHelioviewer movies/
 * ftp://omaftp.oma.be/dist/astro/Verbeeck.F/JHelioviewer states/


Run the notebooks
=================

## Binder
These notebooks (apart from some that require external downloaded data) can be run in the browser using binder.org. 

To launch this click on this link -> [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/SolarOrbiterWorkshop/solo8_tutorials/HEAD)

**Note** It may take a few minutes to load up the first time you launch it.


## Run locally

Here are the steps you'll need to run:
------------------------------------
#### 1. Download the files using git


If you want to run these notebooks locally you can clone this reposity (or fork it and then clone it from your page). To do this run this command:

- ```git clone  https://github.com/SolarOrbiterWorkshop/solo8_tutorials.git```

If you have first forked it then you can run:

- ```git clone  https://github.com/<username>/solo8_tutorials.git```

You can also download these notebooks by clicking on green `code` button on the top right hand side, and then by clicking download zip. 

------------------------------------

#### 2. Create a conda environment

We recommend creating a new conda environment and install the requried packages used in these notebooks. [Here](https://towardsdatascience.com/getting-started-with-python-environments-using-conda-32e9f2779307) is a nice introduction to anaconda environments for those new to the concept, and [here](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf) is a conda cheatsheet which may help too! 

The python packages required to run these workshop notebooks are listed in the `environment.yml` file in this repository. To create a new environment with these packages installed you can open a terminal and type:

- ```conda env create -f environment.yml```

This will then create a new conda environment called `solo8` (this name is listed in the `enviroment.yml` file).

You can then activate this environment by typing:

- ``` conda activate solo8```

Note your prompt should change and now have `solo8` near the start. If you want to list all your conda environments you can type
``` conda info -e```. You should see `base` which is your base enviroment, the `solo8` one, and any others you have created! 


##### 2.2 Updating the environment.yml file
If an update is made to the `environment.yml` file then you will need to type 

- ```conda env update --file environment.yml --prune```

This may be important after you have down a `git pull` (see below 4.)

##### 2.3 Installed new packages in this environment

You can also install new packages in this environment by using `conda install <package>`or by using pip! (`pip install <package<`)


-----------------------------------
### 3. Start a jupyter notebook!

Once you have your environment activated (remember to first type `conda activate solo8`) then in your local `solo8_tutorials` repository type

- ```jupyter notebook ```

This should then open the notebooks in your default browser!

If you are having any issues - just make sure first that you are in the `solo8` environment before you start jupyter.

Happy coding!!

----------------------------------
#### 4. Pulling the most up-to-date version
To make sure that your local repository reflects whats currently here you will need to run a `git pull`. This may be needed is you downloaded the notebooks before the 16th September. This will pull the most up-to-date version of this repository. Before you do this however, you will want to check which remote you have linked to this repository. To find out this you can type this in your local `solo8_tutorials` repository:

- ```git remote -v``` 

and this will list the current remotes
It might looks something like this

```
origin	https://github.com/<username>/solo8_tutorials.git (fetch)
origin	https://github.com/<username>/solo8_tutorials.git (push)
upstream	https://github.com/SolarOrbiterWorkshop/solo8_tutorials.git (fetch)
upstream	https://github.com/SolarOrbiterWorkshop/solo8_tutorials.git (push)
```

what you will want to do is pull the main branch from the one that is linked to `https://github.com/SolarOrbiterWorkshop/solo8_tutorials` - which in this example is `upstream`. Hence to pull the latest version of this repository you would type:

- `git pull upstream main` 

and this will update your local files. 

Otherwise if when you typed `git remote -v` and it looks like
```
origin	https://github.com/SolarOrbiterWorkshop/solo8_tutorials (fetch)
origin	https://github.com/SolarOrbiterWorkshop/solo8_tutorials (push)
```
then you would type

- ``` git pull origin main```

to update your local files. 

