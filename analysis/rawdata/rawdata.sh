#!/bin/bash

# This script automatically generates png files from xml data files.
# These are just initial plots of raw data.
# REQUIREMENTS: assaytools must be installed

for d in $(find ../../data/spectra -mindepth 2 -type d)
do
  #Make figures of xml files in directory $d:
  xml2png --type spectra $d/*.xml
done

mv ../../data/spectra/*/*/*.png spectra

for d in $(find ../../data/singlet -mindepth 0 -type d)
do
  #Make figures of xml files in directory $d:
  xml2png --type singlet_384 $d/*.xml
done

mv ../../data/singlet/*.png singlet
mv ../../data/singlet/*/*.png singlet


