#!/bin/bash

#This script is called by gaussian_force and gaussian_min in order to add the header for a
#Gaussian input file and to run the Gaussian program on that script. The header information
#is added to this script dynamically by executing prepare_gaussian_art.py. The subroutine
#in gaussian_force will set 'force' as a route parameter and the subroutine in gaussian_min
#will set 'opt' as a route parameter.


#*******************************************************************************
# ART Variables
natoms=$1
optimization=$2  #this is 'force' when called from gaussian_force or 'opt' when from gaussian_min
#*******************************************************************************

#*******************************************************************************
# Gaussian input file header and coorLineNumber inserted by prepare_gaussian_art.py
#
#gaussian-header-begin (DO NOT REMOVE) 
header='%mem=8000mb
 %nproc=12
 #rhf/3-21g <PARAM>

name

0 1
'
#Coordinate line where data begins
coorLineNumber=8

#gaussian-header-end
#*******************************************************************************

#*******************************************************************************
# Update header with appropriate options depending on stage in ART

#coorLineNumber=5

updated_header=$()
if [ "$optimization" == "opt" ];
then
    updated_header=$(echo "$header" | sed 's/<PARAM>/nosymm opt/')
elif [ "$optimization" == "force" ];
then
    updated_header=$(echo "$header" | sed 's/<PARAM>/nosymm force/')
fi

# Adds header information to gaussian input file
echo "$updated_header" | cat - art2gaussian.inp > temp && mv temp art2gaussian.inp

#*******************************************************************************


printf "$natoms\n" >>temp.xyz
#TODO title from inp
printf "MOLECULAR TITLE\n" >>temp.xyz

for ((i=1; i<=$natoms; i++)); do
	awk 'FNR=='$coorLineNumber+$i' {print $0}' ./art2gaussian.inp >>temp.xyz
done

# Loads the latest version of Gaussian and calls it through g09
g09 < art2gaussian.inp > test.log

#*******************************************************************************
#Outputs results to a log
printf  "outcoor:\n" >gaussian2art
positionLineNumber=$(sed -n '/Input orientation/=' gaussian.log | tail -1)

for ((i=1; i<=$natoms; i++)); do
	awk 'FNR=='$positionLineNumber+4+$i' {print $0}' ./gaussian.log >>gaussian2art
done

printf  "gaussi: Atomic forces (eV/Ang):\n" >>gaussian2art
forceLineNumber=$(sed -n '/Forces (Hartrees/=' gaussian.log | tail -1)

for ((j=1; j<=$natoms; j++)); do
	awk 'FNR=='$forceLineNumber+2+$j' {print $0}' gaussian.log >>gaussian2art # formats and units are already taken care
done

printf  "gaussi: Final energy (eV):\n" >>gaussian2art

#TODO gets the values of the past ten outputs
for k in {1..10}; do
	grep 'E(' gaussian.log | tail -1 | awk '{printf "gaussi:         Total =  %.10f\n", $5*27.2113838668}' >>gaussian2art
done

#*******************************************************************************
#Outputs results to the gaussian2art file that are read by ART
#
##Number of atoms
#printf  "natoms:\n" >gaussian2art_test
#printf  "$natoms" >gaussian2art_test
#
##Output coordinates
#printf  "outcoor:\n" >gaussian2art_test
#positionLineNumber=$(sed -n '/Input orientation/=' gaussian.log | tail -1)
#
##Energy
#printf  "energy:\n" >>gaussian2art_test


#Deprecated

# for b3lyp
#sed -i "s/rhf/b3lyp/g" art2gaussian.inp # change method
#sed -i "s/3-21g/6-311++G(d,p)/g" art2gaussian.inp # change basis set

# Maps from numerical value to letter for atoms
# sed 's/\([^ ]*\) [0-9]*[1-9][0-9]* /\1 C /' temp.xyz > output.xyz
# sed 's/\([^ ]*\) [6]*[6][6]* /\1 C /; s/\([^ ]*\) [1]*[1][1]* /\1 H /' temp.xyz > ethane_output.xyz
# for temp.xyz STARTS