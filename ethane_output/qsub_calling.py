#!/usr/bin/python
#PBS -l nodes=4:ppn=12
#PBS -l walltime=01:00:00
#PBS -A xyz-123-ab
#PBS -j oe
#PBS -N ethane_output

import os
import subprocess

os.chdir(os.getenv('PBS_O_WORKDIR', '/global/scratch/bhu2108'))

os.environ['IPATH_NO_CPUAFFINITY'] = '1'

os.environ['OMP_NUM_THREADS'] = '12'

subprocess.call("qsubcmd", shell=True)
