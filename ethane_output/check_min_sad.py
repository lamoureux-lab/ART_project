import argparse
parser = argparse.ArgumentParser(description = 'Create an input and submission file')
parser.add_argument('filename', type = str, help = 'Name of the input file')
parser.add_argument('logfile', type = str, help = 'Name of the output file')
args = parser.parse_args()

def create_input_file(filename):
	with open(filename) as f, open(filename + '.com','w') as m:
		for i in range(0,3):
			_= f.readline()
		if filename.startswith("mi"):
			m.write("%" + "chk" + " = " + filename + ".chk\n# opt freq hf/6-31g\n\nEthane Optimization\n\n0 1\n")
		elif filename.startswith("sad"):
			m.write("%" + "chk" + " = " + filename + ".chk\n# opt(qst3) freq hf/6-31g\n\nEthane Optimization\n\n0 1\n")
		for line in f:
			m.write(line)
		m.write("\n\n\n\n\n\n\n\n")
		m.close()

def create_submission_file(filename):
	with open(filename + '.sub','w') as n:
		n.write('''#!/bin/bash \n #PBS -S /bin/bash \n #PBS -l nodes=1:ppn=4 \n #PBS -l mem=1800MB \n #PBS -l walltime=5:00:00 \n #PBS -N Ethane_Opt \n # Adjust the mem and ppn above to match the requirements of your job \n # Sample Gaussian job script \n cd $PBS_O_WORKDIR \n echo "Current working directory is `pwd`" \n echo "Running on `hostname`" \n echo "Starting run at: `date`" \n # Set up the Gaussian environment using the module command: \n module load gaussian \n # Run Submission \n g09 ''' + filename + '''.com \n \n \n''')
		n.close()

def check_min_or_sad(logfile):
	with open(logfile)as f:
	    frequency = []
	    for line in f:
	        if line.startswith(" Frequencies"):
	            frequency.append(line)
	            print(line)

	for line in frequency:
		freq = line.split()
		check = float(freq[2])
		if check < 0:
			print("Saddle")

if __name__ == '__main__':
	create_input_file(args.filename)
	create_submission_file(args.filename)
	check_min_or_sad(args.logfile)
