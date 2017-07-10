
## ART Setup Process    
 
1. Apply for a Compute Canada Account:
https://ccdb.computecanada.ca
 
2. Apply for a Westgrid Consortium account:
	Following approval of the Compute Canada account, login and apply for a WestGrid 
  Consortium account for access to GREX
 
3. Apply for Gaussian access:
https://www.westgrid.ca/support/software/gaussian/license_terms
 
 
## GREX Basics 
 
Access GREX server:
    `ssh <westgrid_username>@grex.westgrid.ca`
 
This will log you into your home directory
	`/home/<username>/`
 
Running and testing code:
	Move files to 
    `/global/scratch/<username>/`
 
## Running 
 
Ensure that the ART_project directory is in `/global/scratch/<username>/` and not `/home/<username>/` as this is a heavy process. 

Once there, compile the ART code by running 
    `make /global/scratch/<username>/ART_project/gaussian_art_v8 artdft`
 
In the sample file directory, execute the following:
	`csh siestart.sh > output.log`

#TODO 
What is this for/ does it belong in the readme?  
grep 'E(' test.log | tail -1 | awk '{printf "siesta:         Total =  %.10f\n", $5*27.2113838668}' >>log

Change the following configuration values in accordance to the sample:
1: execute_gaussian.sh => number of atoms natoms
2: change siestart.sh => NATOMS, TYPE OF ATOMS and other as require 
 
 
## Troubleshooting
 
Server status: 
https://www.westgrid.ca/support/systems/Grex/system_status
 
For additional questions, contact support@westgrid.ca

##File structure 
clear.sh is a simple script that remove all generated files an resets 

