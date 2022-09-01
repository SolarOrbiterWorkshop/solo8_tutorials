# SolO-HI tutorial

This holds the SolO-HI tutorial.

Sample data for this tutorial can be accessed at:

https://app.box.com/s/y7s5fdgof6a6p1xj6jusy56y3zzwthn2

All routines needed to combine files from individual SoloHI tiles into a mosaic and generate images, movies, and J-maps(elongation vs. time plots) are given in shifits2grid.py

The PATH variable in line 2 of the Jupyter notebook should be set to the where the data is located. The data should ideally be stored in separate subdirectories by day, but can also be combined into a single folder. If this is done however, the DATES variable should be updated to reflect the data subdirectory

To get the movie generation tools to function, the user must have FFMPEG installed and:

plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'

in line 9 must be updated to reflect the path to your local FFMPEG binary. Different writers (i.e. avconv) can also be used by modifying the FuncAnimation save line in shifits2grid.py (line 402)

These routines have primarily been developed to work on Level 2 and higher data products, and this is the example presented in the notebook. If you want to use L1 data for whatever reason, the code has been provided that will provide a (greatly simplified) calibration process to remove the bias, correct detector linearity and normalize for exposure time. The 8 IDL .sav files provided give needed data to perform this process. The SoloHI team encourages the use of Level 3 or higher products.

If you have any questions, or need help accessing the data, contact Phillip Hess (phillip.hess@nrl.navy.mil) for help.
