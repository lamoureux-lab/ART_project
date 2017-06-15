#!/bin/bash

#-------------------------------------------------------------------------------
# ART Variables
natoms=$1
optimization=$2  #this is 'force' when generated from gaussian_force or 'opt' when generated from gaussian_min
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Gaussian input file header and coorLineNumber inserted by run_gaussian_artdft.py
#
#gaussian-header-begin (DO NOT REMOVE) 
header='%mem=8000mb
 %nproc=12
 #rhf/3-21g <OPTION>

name

0 1
'
#Coordinate line where data begins
coorLineNumber=8

#gaussian-header-end
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Update header with appropriate options depending on stage in ART

updated_header=$()
if [ "$optimization" == "opt" ];
then
    updated_header=$(echo "$header" | sed 's/<OPTION>/nosymm opt/')
elif [ "$optimization" == "force" ];
then
    updated_header=$(echo "$header" | sed 's/<OPTION>/nosymm force/')
fi

# Adds header information to gaussian input file
echo "$updated_header" | cat - art2gaussian.inp > temp && mv temp art2gaussian.inp

#-------------------------------------------------------------------------------


printf "$natoms\n" >>temp.xyz
printf "MOLECULAR TITLE\n" >>temp.xyz


for ((i=1; i<=$natoms; i++)); do
	awk 'FNR=='$coorLineNumber+$i' {print $0}' ./art2gaussian.inp >>temp.xyz
done

# Loads the latest version of Gaussian and calls it through g09
g09 < art2gaussian.inp > gaussian.log

#----------------------------------------------------
#Outputs results to a log
printf  "outcoor:\n" >log
positionLineNumber=$(sed -n '/Input orientation/=' gaussian.log | tail -1)

for ((i=1; i<=$natoms; i++)); do
	awk 'FNR=='$positionLineNumber+4+$i' {print $0}' ./gaussian.log >>log
done

printf  "gaussi: Atomic forces (eV/Ang):\n" >>log
forceLineNumber=$(sed -n '/Forces (Hartrees/=' gaussian.log | tail -1)

for ((j=1; j<=$natoms; j++)); do
	awk 'FNR=='$forceLineNumber+2+$j' {print $0}' gaussian.log >>log # formats and units are already taken care
done

printf  "gaussi: Final energy (eV):\n" >>log

for k in {1..10}; do
	grep 'E(' gaussian.log | tail -1 | awk '{printf "gaussi:         Total =  %.10f\n", $5*27.2113838668}' >>log
done
#----------------------------------------------------

#Deprecated

# for b3lyp
#sed -i "s/rhf/b3lyp/g" art2gaussian.inp # change method
#sed -i "s/3-21g/6-311++G(d,p)/g" art2gaussian.inp # change basis set

# Maps from numerical value to letter for atoms
# sed 's/\([^ ]*\) [0-9]*[1-9][0-9]* /\1 C /' temp.xyz > output.xyz
# sed 's/\([^ ]*\) [6]*[6][6]* /\1 C /; s/\([^ ]*\) [1]*[1][1]* /\1 H /' temp.xyz > ethane_output.xyz
# for temp.xyz STARTS