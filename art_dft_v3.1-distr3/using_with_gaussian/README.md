
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
 
Ensure that the ART_project directory is in `/global/scratch/<username>/` and not `/home/<username>/` as this is 
a heavy process. 

Once there, compile the ART program by navigating to the ART source directory and running 
    `make artdft`
	
To run the simulations, navigate to a directory containing Gaussian input files and run:
    `python <path_to_scripts>/run_gaussian_art.py`

Review program usage:
    `python <path_to_scripts>/run_gaussian_art.py -h`

By default, this will load all `<structure>.inp` Gaussian input files in the `/work` directory, and prepare
output directories `/work/<structure1>`, `/output/<structure2>`, etc. with a copy of the execution
scripts to run ART with Gaussian. A submission script from each output directory is then submitted using -qsub to 
run these jobs in parallel. All min and saddle results are stored in the output directories. 

Review your queued jobs and their progress:
    `qstat -u <username>`

For other instruction on how to use qsub:
https://wikis.nyu.edu/display/NYUHPC/Tutorial+-+Submitting+a+job+using+qsub#Tutorial-Submittingajobusingqsub-what_is_qsub
http://beige.ucs.indiana.edu/I590/node35.html


 
## Troubleshooting
 
Server status: 
https://www.westgrid.ca/support/systems/Grex/system_status
 
For additional questions, contact support@westgrid.ca

