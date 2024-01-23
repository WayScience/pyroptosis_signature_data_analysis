# Generating publication ready figures from the data

## Requirements
- FIJI (ImageJ) version 1.54f
    - Download [FIJI](https://imagej.net/software/fiji/downloads)
FIJI stands for FIJI is just ImageJ.
FIJI is a distribution of ImageJ which includes a lot of plugins and macros.

## Before you start
Before running the scripts, to generate figure5, we ran FIJI macros to adjust the lookup tables of each image.
We change red to magenta, keep green as green, and change blue to cyan.
These changes make the images for accessible to people with color blindness.

The macro can be found under the folder `figures/5.figure5/imageJ_macros/`.

Please run the macros before running the scripts that assemble the figures.

Please note that we used version 1.54f of FIJI.

Open FIJI and drag and drop the macro file into the FIJI window.
Then click on `Run` and the macro will run.
The macro is written for user interface and will ask you to specify the input and output folders.
