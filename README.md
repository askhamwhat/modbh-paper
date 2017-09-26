# modbh-paper

The routines in this directory were
used to make the figures in the paper
"A stabilized separation of variables 
method for the modified biharmonic
equations"

## About

This repository is intended to be a
snapshot of the code as it was used
in the paper. A separate repository
will be created for housing the
underlying code and any further updates

## How to use

To generate a makefile, run the script
configure_makefile.sh. It should then
be possible to simply type make to
generate the figures. The figures will end
up in the folder named "fig".

## Requirements

To generate the figures you must have Python
and MatPlotLib. As of yet, this code has
only been tested on a Linux platform with
the gfortran compiler.