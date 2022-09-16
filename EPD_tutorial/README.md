# EPD tutorial

Example notebook for the Python data loader [solo-epd-loader](https://github.com/jgieseler/solo-epd-loader) for Solar Orbiter's (SolO) Energetic Particle Detector (EPD), provided through the [SERPENTINE project](https://serpentine-h2020.eu).

## Run in Binder
This notebook can be run online in Binder following this link [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/SolarOrbiterWorkshop/solo8_tutorials/HEAD/?filepath=EPD_tutorial/epd_data_loader.ipynb)

## Installation


[solo-epd-loader](https://github.com/jgieseler/solo-epd-loader) requires Python >= 3.6.

It can be installed either from [PyPI](https://pypi.org/project/solo-epd-loader/) using:

``` bash
pip install solo-epd-loader
```   

or from [conda](https://anaconda.org/conda-forge/solo-epd-loader/) using:

``` bash
conda install -c conda-forge solo-epd-loader
```

In addition, to run the Jupyter Notebook, you need the [jupyter package installed](https://jupyter.org/install).


## Run 
1. Open your terminal/command line/Anaconda prompt.
2. In the terminal, navigate to the folder that contains the file `epd_data_loader.ipynb`
3. Run `jupyter notebook epd_data_loader.ipynb`
4. Your standard web-browser should now open the Jupyter Notebook.


## Additional tools
Within the [SERPENTINE](https://serpentine-h2020.eu) project, additional tools are developed that can be found at https://github.com/serpentine-h2020/serpentine/tree/main/notebooks/sep_analysis_tools. These Jupyter Notebooks include for example functions to generate dynamic spectra or determine Solar Energetic Particle (SEP) onset times utilizing EPD observations.


## Acknowledgements
This project has received funding from the European Union’s Horizon 2020 research and innovation programme under grant agreement No 101004159.
