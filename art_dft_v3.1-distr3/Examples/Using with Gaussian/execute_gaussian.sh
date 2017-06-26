#!/bin/bash

#This script is called by gaussian_force and gaussian_min in order to add the header for a
#Gaussian input file and to run the Gaussian program on that script. The header information
#is added to this script dynamically by executing prepare_gaussian_art.py. The subroutine
#in gaussian_force will set 'force' as a route parameter and the subroutine in gaussian_min
#will set 'opt' as a route parameter.
gaussian2art='gaussian2art'
art2gaussian='art2gaussian.inp'
gaussian_output='gaussian.log'

#*******************************************************************************
# ART Variables
natoms=$1
optimization=$2  #this is 'force' when called from gaussian_force or 'opt' when from gaussian_min

#*******************************************************************************
#*******************************************************************************

# Gaussian input file header data, including the coorLineNumber and title information
# are inserted dynamically by prepare_gaussian_art.py between the gaussian-header-begin
# and gaussian-header-end flags.
#
#gaussian-header-begin (DO NOT REMOVE)

#gaussian-header-end (DO NOT REMOVE)

#*******************************************************************************
#*******************************************************************************
# Update header with appropriate options depending on stage in ART

updated_header=$()
if [ "$optimization" == "opt" ];
then
    updated_header=$(echo "$header" | sed 's/<PARAM>/nosymm opt/')
elif [ "$optimization" == "force" ];
then
    updated_header=$(echo "$header" | sed 's/<PARAM>/nosymm force/')
fi

# Adds header information to gaussian input file
echo "$updated_header" | cat - "$art2gaussian" > temp && mv temp "$art2gaussian"

#*******************************************************************************
#*******************************************************************************

# Updates running record of atom coordinates throughout program execution

printf "$natoms\n" >>temp.xyz
printf "$title\n" >>temp.xyz

for ((i=1; i<=$natoms; i++)); do
	awk 'FNR=='$coorLineNumber+$i' {print $0}' "$art2gaussian" >>temp.xyz
done


#*******************************************************************************
#*******************************************************************************
# Loads the latest version of Gaussian and calls it through g09, saving the output
# to gaussian.log
g09 < $art2gaussian > $gaussian_output

#*******************************************************************************
#*******************************************************************************
# Prepares the gaussian results for ART in the gaussian2art file format

#Number of atoms
printf  "natoms:\n" >$gaussian2art
printf  "$natoms\n" >>$gaussian2art

# Coordinates
printf  "outcoor:\n" >>$gaussian2art
positionLineNumber=$(sed -n '/Input orientation/=' gaussian.log | tail -1)

for ((i=1; i<=$natoms; i++)); do
	TEMP=$(awk 'FNR=='$positionLineNumber+4+$i' {print $0}' gaussian.log)
	TEMP=$(echo $TEMP | sed 's/[\S]+//')
	echo "$TEMP" | cut -d ' ' -f4- >>$gaussian2art
done

# Energy
printf  "energy:\n" >>$gaussian2art
grep 'E(' "$gaussian_output" | tail -1 | awk '{printf "%.10f\n", $5*27.2113838668}' >>$gaussian2art

# Forces
forceLineNumber=$(sed -n '/Forces (Hartrees/=' gaussian.log | tail -1)
# Prepares force information when available
if [ "$forceLineNumber" != "" ]
then
    printf  "forces:\n" >>$gaussian2art
    for ((j=1; j<=$natoms; j++)); do
        TEMP=$(awk 'FNR=='$forceLineNumber+2+$j' {print $0}' gaussian.log)
        TEMP=$(echo $TEMP | sed 's/[\S]+//')
        echo "$TEMP" | cut -d ' ' -f3- >>$gaussian2art
    done
fi


#*******************************************************************************

##*******************************************************************************
##Outputs results to the gaussian2art file that are read by ART (OLD STRUCTURE)
#printf  "outcoor:\n" >gaussian2art
#positionLineNumber=$(sed -n '/Input orientation/=' gaussian.log | tail -1)
#
#for ((i=1; i<=$natoms; i++)); do
#	awk 'FNR=='$positionLineNumber+4+$i' {print $0}' ./gaussian.log >>gaussian2art
#done
#
#printf  "gaussi: Atomic forces (eV/Ang):\n" >>gaussian2art
#forceLineNumber=$(sed -n '/Forces (Hartrees/=' gaussian.log | tail -1)
#
#for ((j=1; j<=$natoms; j++)); do
#	awk 'FNR=='$forceLineNumber+2+$j' {print $0}' gaussian.log >>gaussian2art # formats and units are already taken care
#done
#
#printf  "gaussi: Final energy (eV):\n" >>gaussian2art
#

#for k in {1..10}; do
#	grep 'E(' gaussian.log | tail -1 | awk '{printf "gaussi:         Total =  %.10f\n", $5*27.2113838668}' >>gaussian2art
#done


#Deprecated

# for b3lyp
#sed -i "s/rhf/b3lyp/g" art2gaussian # change method
#sed -i "s/3-21g/6-311++G(d,p)/g" art2gaussian # change basis set

# Maps from numerical value to letter for atoms
# sed 's/\([^ ]*\) [0-9]*[1-9][0-9]* /\1 C /' temp.xyz > output.xyz
# sed 's/\([^ ]*\) [6]*[6][6]* /\1 C /; s/\([^ ]*\) [1]*[1][1]* /\1 H /' temp.xyz > ethane_output.xyz
# for temp.xyz STARTS