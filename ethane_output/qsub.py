import argparse
parser = argparse.ArgumentParser(description = 'Create a submission file')
parser.add_argument('filename', type = str, help = 'Name of the file')
args = parser.parse_args()

def create_submission_file(filename):
	with open(filename + '.sub','w') as f:
		f.write('''#!/bin/bash \n #PBS -S /bin/bash \n #PBS -l nodes=1:ppn=4 \n #PBS -l mem=1800MB \n #PBS -l walltime=5:00:00 \n #PBS -N Ethane_Opt \n # Adjust the mem and ppn above to match the requirements of your job \n # Sample Gaussian job script \n cd $PBS_O_WORKDIR \n echo "Current working directory is `pwd`" \n echo "Running on `hostname`" \n echo "Starting run at: `date`" \n # Set up the Gaussian environment using the module command: \n module load gaussian \n # Run Submission \n g09 ''' + filename + '''.com''')

if __name__ == '__main__':
	create_submission_file(args.filename)

