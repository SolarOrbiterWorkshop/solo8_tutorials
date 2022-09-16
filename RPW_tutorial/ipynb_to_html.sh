#!/bin/bash

# Create output directory (empty)
rm -rf ./html
mkdir -p ./html

# Convert Jupyter notebboks to html files
for current_file in `ls *.ipynb`;do
  jupyter nbconvert --to html $current_file
done

# Move html files in html/ folder
mv *.html html/.

