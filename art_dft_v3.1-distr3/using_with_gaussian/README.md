
## GREX Setup Process    
 
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
 
For quick tests to ensure functionality, execute the following in the sample file directory:
	`csh gaussian_art.sh > output.log`
	
For actual simulations, submit the job to the GREX queue using:
    `qsub -N art-test gauss_grex.sub`

Review your queued jobs and their progress:
    `qstat -u <username>`

For other instruction on how to use qsub:
https://wikis.nyu.edu/display/NYUHPC/Tutorial+-+Submitting+a+job+using+qsub#Tutorial-Submittingajobusingqsub-what_is_qsub
http://beige.ucs.indiana.edu/I590/node35.html

Change the following configuration values in accordance to the sample:
1: execute_gaussian.sh => number of atoms natoms
2: change gaussian_art.sh => NATOMS, TYPE OF ATOMS and other as require 
 
 
## Troubleshooting
 
Server status: 
https://www.westgrid.ca/support/systems/Grex/system_status
 
For additional questions, contact support@westgrid.ca

