#!/bin/bash

# Create output directory (empty)
rm -rf html
mkdir html

# Convert Jupyter notebboks to html files
jupyter nbconvert --to html 1_about_rpw_data.ipynb
jupyter nbconvert --to html 2_get_rpw_data.ipynb
jupyter nbconvert --to html 3a_rpw_tnr_data_tutorial.ipynb
jupyter nbconvert --to html 3b_rpw_tds_data_tutorial.ipynb

# Move html files in html/ folder
mv *.html html/.

