#!/bin/bash

natoms=8

sed -i '1d' art2gaussian.inp # removes chk point line
#sed -i "s/0 1/1 1/g" art2gaussian.inp # change charge and multiplicity

# for b3lyp
#sed -i "s/rhf/b3lyp/g" art2gaussian.inp # change method
#sed -i "s/3-21g/6-311++G(d,p)/g" art2gaussian.inp # change basis set

# Maps from numerical value to letter for atoms
# sed 's/\([^ ]*\) [0-9]*[1-9][0-9]* /\1 C /' temp.xyz > output.xyz
# sed 's/\([^ ]*\) [6]*[6][6]* /\1 C /; s/\([^ ]*\) [1]*[1][1]* /\1 H /' temp.xyz > ethane_output.xyz
# for temp.xyz STARTS

printf "$natoms\n" >>temp.xyz
printf "MOLECULAR TITLE\n" >>temp.xyz

#Coordinate line where the data begins
coorLineNumber=5

for ((i=1; i<=$natoms; i++)); do
	awk 'FNR=='$coorLineNumber+$i' {print $0}' ./art2gaussian.inp >>temp.xyz
done

# for temp.xyz ENDS

sed -i '1 i %nproc=12' art2gaussian.inp 
sed -i '1 i %mem=8000MB' art2gaussian.inp

# Loads the latest version of Gaussian and calls it through g09
g09 < art2gaussian.inp > test.log
#more test.log >>log1

printf  "outcoor:\n" >log
positionLineNumber=$(sed -n '/Input orientation/=' test.log | tail -1)

for ((i=1; i<=$natoms; i++)); do
	awk 'FNR=='$positionLineNumber+4+$i' {print $0}' ./test.log >>log
done 

printf  "siesta: Atomic forces (eV/Ang):\n" >>log
forceLineNumber=$(sed -n '/Forces (Hartrees/=' test.log | tail -1)

for ((j=1; j<=$natoms; j++)); do
	awk 'FNR=='$forceLineNumber+2+$j' {print $0}' test.log >>log # formats and units are already taken care
done

printf  "siesta: Final energy (eV):\n" >>log

for k in {1..10}; do
	grep 'E(' test.log | tail -1 | awk '{printf "siesta:         Total =  %.10f\n", $5*27.2113838668}' >>log
done

