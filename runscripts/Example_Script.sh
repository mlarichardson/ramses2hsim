#!/bin/bash

# CODE: ramses2hsim
# AUTHOR: Mark Richardson (Mark.Richardson@McDonaldInstitute.ca)
# DATE: 2020-08-20
# PURPOSE: This pipeline takes a RAMSES (Teyssier 2002) output, extracts a region, and converts it to
#          a data cube (arcsec x arcsec x microns) that can be used as an input for HSIM (Zieleniewski et al. 2016).
# PAPER: Please see Richardson, Routledge, Thatte, et al. 2020, MNRAS. 



# PREREQUISITS:
#  - Have an output_0000N file from ramses on which you want to run this pipeline
#  - Have located a galaxy with position XC YC ZC and radius RAD in grid units
#
#  - Have software:
#     - amr2Intcube
#     - intcube2velcube
#     - convert_to_fits.py
#
#     - Makefile
#
#     - dependency: ramses_info.f90, binaryio.py, numpy, astropy, scipy
#
#  - Assuming here you are running pipeline in a subdirectory (perhaps titled PIPELINE), such that the
#    path to the output info file is: ../output_0000N/info_0000N.txt
#

# STEP 0: COMPILE CODE

make all




# STEP 1: CONVERT RAMSES OUTPUT REGION TO UNIFORM RESOLUTION INTERMEDIATE CUBE

# Run amr2Intcube in desired (or all three) projection directions:
#  In this example, N = 152 (-inp = input name: -inp ../output_00152)
#                   -xc XC = 0.70517813
#                   -yc YC = 0.33494881
#                   -zc ZC = 0.34129348
#                   -rc RAD = 1.2253e-3 (this sim is at z=3 with boxsize=12.5cMpc --> R = 3.83 pkpc.
#
#  Extra arguments:
#    -out >> Name of Output file. I keep .bin to indicate it's a binary output.
#    -dir >> Direction to project along. Default is z
#    -typ >> method to determine Ha. Choices are sfr_Ha and Ost_Ha.
#                  Other quantities can be plotted. All listed in ramses_info emissivity subroutine.
#    -lma >> maximum level to project to. Note, in this example simulation, at z=3, max level reached is 18, dx = 10pc. -lma means only save data to 10pc.
#
#  Other possible arguments:
#
#
#  Result:
#   Constructs a 3d uniform mesh table, where each cell has variables:
#     -Ha emission
#     -Metallicity
#     -LOS velocity
#     -LOS velocity dispersion
#
#   **NB: Sometimes depending on RAD and LMAX, this may require significant memory.
#         You may get around this by setting: ulimit -s unlimited.
#
XC=0.70517813
YC=0.33494881
ZC=0.34129348
RAD=1.2253e-3
NUM=00152
trimNum="`echo $NUM | sed 's/^0*//'`"
DIR=x
LMAX=17
METH="sfr_Ha"
amr2Intcube -inp ../output_$NUM -out IntCube_${DIR}_${trimNum}_l${LMAX}_${METH}.bin -xc $XC -yc $YC -zc $ZC -rc $RAD  -dir $DIR -typ "$METH" -lma $LMAX


STEP 2: CONVERT UNIFORM RESOLUTION INTERMEDIATE CUBE TO VELOCITY CUBE

# Run intcube2velcube
#  In this example:
#    -inp INPUT_NAME = output of amr2Intcube = IntCube_x_152_l17_sfr_Ha.bin
#    -out OUTPUT_NAME
#    -fd FDUST, = dust_fraction, the fraction of metals, by mass, residing in dust.
#
#  Other Parameters:
#   -vmin, vmax, nv
#
#  Result:
#    Constructs a 2 spatial dimension x 1 velocity dimension cube, where each cell has the Ha luminosity.
#
FDUST=0.01
intcube2velcube -inp IntCube_${DIR}_${trimNum}_l${LMAX}_${METH}.bin -out Velcube_${DIR}_${trimNum}_l${LMAX}_${METH}_fd${FDUST}.bin -fd $FDUST





# STEP 3: CONVERT VELOCITY CUBE TO WAVELENGTH SPACE

# Run convert_to_fits.py
#  This code takes the rest-frame cube from amr2velcube_usingIntCube, and redshifts it
#  to a desired redshift, converting the velocity dimension into wavelength dimension centred on a given wavelength (redshifted).
#  Currently this needs to be editted by hand before calling:
#
#  N.B: You need to set filenames to the list of basenames (so not including .bin) of output files from
#  amr2velcube_usingIntCube.
#
# You can get this like:
# >> ls -1 Velcube_* | sed 's/^/"/' | sed 's/.bin/",/'

#  You also need to set the scale factors, aexps, to the values you want to redshift the vel cube into.
#

python convert_to_fits.py
