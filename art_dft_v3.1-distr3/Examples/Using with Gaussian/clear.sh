#!/bin/bash
#This script can be used to clear files produced after a run
#TODO change output files to a distinct folder

sed -i "/Counter:/c\Counter:      1000" filecounter
rm art2siesta CLOCK fdf-*.log INPUT_TMP.* log.file.* output.log siesta2art Si+I.* test_*
rm sad1*
rm art2gaussian gaussian2art.log temp.chk log events.list test.log temp.xyz
