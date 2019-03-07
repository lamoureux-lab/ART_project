import argparse
from os import getcwd

parser = argparse.ArgumentParser(description = 'ART coupled with CP2K')

parser.add_argument('-k', '--keyword', help = 'input file')

args =  parser.parse_args()

cp2k_input_file = 'cp2k_input.inp'

with open ('../ethane2.inp') as f:
    with open (cp2k_input_file, 'w+') as m:
        for line in f:
            if 'RUN_TYPE' in line:
                if args.keyword == 'opt':
                    line = line.replace(line.split()[1], 'GEO_OPT')
                elif args.keyword == 'ef':
                    line = line.replace(line.split()[1], 'ENERGY_FORCE')
            m.write(line)
            
        

