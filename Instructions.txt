
 ## Running  

Once in the ART source directory, compile ART by running 
    `make artdft`
	
To run the simulations, navigate to a directory containing Gaussian input files and run:
    `python <path_to_scripts>/artgau.py`

Review program usage:
    `python <path_to_scripts>/artgau.py -h`

By default, this will load all `<structure>.inp` Gaussian input files in the `/work` directory, and prepare
output directories `/work/<structure1>`, `/output/<structure2>`, etc. with a copy of the execution
scripts to run ART with Gaussian. A submission script from each output directory is then submitted using "sbatch' to 
run these jobs in parallel. All min and saddle results are stored in the output directories. 

Review your queued jobs and their progress:
    `squeue -u <username>`
 

