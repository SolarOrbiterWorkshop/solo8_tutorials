# EUI tutorial

The EUI tutorial is spread out over (currently) four different Jupyter Notebooks. It's highly suggested to start at the beginning, with the notebook labeled 1_introduction. The notebooks are not stand-alone, but make use of an EUI dataset that can be downloaded as a separate zip file from the ROB cloud system [here](https://cloud-as.oma.be/index.php/s/jzZ5qHPqzb8zMrN). 

This data.zip file (273 MiB) should be extracted in the same directory as the Jupyter Notebooks, creating several 'data_' directories next to the .ipynb notebook files. When uncompressed, the zip file takes up approximately 1.5 GiB, and contains two SQLite databases as well as several SolO EUI and other FITS files that are used in the Jupyter Notebooks.

In order to run the Notebooks, the following Python packages need to be installed. If you're running anaconda, then from a conda prompt, type "pip install ..." replacing ... with the names below :
- matplotlib
- sunpy
- astropy
- opencv-python
- pandas
- sqlite3
- scipy
- soloEUI ( pip install soloEUI --extra-index-url https://gitlab-as.oma.be/api/v4/projects/581/packages/pypi/simple )

If you don't already have Jupyter Notebook installed, then from  the same conda prompt, you can install Jupyter Notebook by typing "pip install notebook". To start Jupyter Notebook, simply type "jupyter notebook", which will open a browser you can use to open Jupyter Notebook files.

Please note that some of the tutorials are quite memory intensive, in particular the 4_euv_mosaic notebook may fail to run completely on systems with less than 16 GiB of memory.

