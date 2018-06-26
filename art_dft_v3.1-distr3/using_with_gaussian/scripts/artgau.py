from os.path import dirname, join, relpath, isfile, exists
from os import makedirs
import argparse
import re
from subprocess import call

parser = argparse.ArgumentParser(description = 'ART coupled with Gaussian')

parser.add_argument('-f', '--gaussian_input', nargs = '+', help = 'Gaussian input file')
parser.add_argument('-m', '--mem', help = 'memory allocation, e.g. 16G', default='16G')
parser.add_argument('-n', '--nodes', help = 'processor nodes, e.g. 2', default = '2')
parser.add_argument('-t', '--walltime', help = 'walltime, e.g. 20:00:00', default = '20:00:00')
parser.add_argument('-g', '--gpus', help = 'number of gpus required, e.g. 2', default ='2')
parser.add_argument('-c', '--cpus', help = 'number of cpus per task, e.g. 16', default ='16')
parser.add_argument('-ne', '--max_events', help = 'number of ART events, e.g. 1', default ='1')
parser.add_argument('-type', '--type_events', help = 'type of ART event, e.g. global or local', default ='global')
parser.add_argument('-r', '--radius', help = 'radius of initial deformation, e.g. 3.0', default ='3.0')
parser.add_argument('-cat', '--centre_atm', help = 'central atom for local event, e.g. 1', default ='1')
parser.add_argument('-eigen', '--eigen_thresh', help = 'eigenvalue threshold, e.g. -0.2', default ='-0.2')
parser.add_argument('-force', '--force_thresh', help = 'exit force thresh, e.g. 0.1', default ='0.1')
parser.add_argument('-inc', '--inc_size', help = 'increment size', default ='0.1')
parser.add_argument('-init', '--init_step', help = 'initial step size, e.g. 0.1', default ='0.1')
parser.add_argument('-sym', '--sym_break', help = 'symmetry break distance, e.g. 0.2', default ='0.2')
parser.add_argument('-pb', '--max_perp_basin', help = 'max perp steps in basin, e.g. 3', default ='3')
parser.add_argument('-k', '--min_k', help = 'minimum number of k-steps, e.g. 3', default ='3')
parser.add_argument('-bf', '--basin_factor', help = 'factor multiplying increment size for leaving the basin, e.g. 2.1', default ='2.1')
parser.add_argument('-l', '--lanczos', help = 'number of iterations in the Lanczos loop, e.g. 5', default ='5')
parser.add_argument('-am', '--activ_maxiter', help = 'max iterations during activation, e.g. 400', default ='400')
parser.add_argument('-d', '--delta_thresh', help = 'energy threshold during Lanczos e.g. 4.0', default ='4.0')
parser.add_argument('-pa', '--max_perp_activ', help = 'max perp steps during activation e.g 5', default ='5')
parser.add_argument('-fp', '--force_thresh_perp', help = 'force thresh perp relax, e.g. 0.5', default ='0.5')
parser.add_argument('-bm', '--basin_maxiter', help = 'max iterations in basin, e.g. 20', default ='20')
parser.add_argument('-xyz', '--write_xyz', help = 'write or not xyz files -- true or false', default ='true')
parser.add_argument('-search', '--search_strat', help = 'search strategy, e.g. 0, 1 or 2', default ='0')

args =                        parser.parse_args()

submission_script_template =  join(dirname(relpath(__file__)), 'gauss_graham.sub')
exec_script_template =        join(dirname(relpath(__file__)), 'execute_gaussian.py')
create_header_template =      join(dirname(relpath(__file__)), 'create_header.py')
gauss_art_template =          join(dirname(relpath(__file__)), 'gaussian_art.sh')

submission_script =          'gauss_graham.sub'

exec_script =                'execute_gaussian.py'

create_header =              'create_header.py'

gauss_art =                  'gaussian_art.sh'

refconfig =                  'refconfig.dat'

filecounter =                'filecounter'

files_to_run =                args.gaussian_input

my_dict =                     {}

for each_file in files_to_run:
    my_dict[each_file] = re.split(r'[.]', each_file)[0]

def create_directory():

    for each_file, directory in my_dict.iteritems():

        try:
            if not exists(directory):
                makedirs(directory)
        except OSError:
            print ('Error: Creating directory. ' +  directory)


def create_submission_file():

    for directory in my_dict.itervalues():        
        with open(submission_script_template) as f:
            with open(join(directory,submission_script), 'w+') as m:
                for line in f:
                    if 'mem' in line:
                        line = line.replace(line.split('=')[1], args.mem)
                    if 'nodes' in line:
                        line = line.replace(line.split('=')[1], args.nodes)
                    if 'time' in line:
                        line = line.replace(line.split('=')[1], args.walltime)
                    if 'gres' in line:
                        line = line.replace(line.split('=')[1], 'gpu:' + args.gpus)
                    if 'cpus' in line:
                        line = line.replace(line.split('=')[1], args.cpus)

                    m.write(line + '\n')
               
def create_exec_script():
 
    for directory in my_dict.itervalues():
        with open(exec_script_template) as f:
            with open(join(directory, exec_script), 'w+') as m:
                for line in f:
                    m.write(line)

def create_header_script():

    for each_file, directory in my_dict.iteritems():
        with open(create_header_template) as f:
            with open(join(directory,create_header),'w+') as m:
                    for line in f:
                        if "input_file =" in line:
                            line = line.replace("input_file =", "input_file = " + "'" + "../" + each_file + "'")
                        m.write(line)

def calculate_atoms():
    atom_list = []
    for each_file in files_to_run:
        natoms = 0
        with open(each_file) as f:
            for line in f:
                if (re.findall(r'\s\d\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)):
                    natoms = natoms + 1
            
            atom_list.append(natoms)

    zipped_atomic_list = zip(files_to_run, atom_list)

    return zipped_atomic_list


def create_gauss_art(zipped_atomic_list):
 
    for each_file, natoms in zipped_atomic_list:
        with open(gauss_art_template) as f:
            with open(join(each_file.split('.')[0], gauss_art), 'w+') as m:
                for line in f:                     
                    if 'NATOMS' in line:
                        line = line.replace(line.split('#')[0].split()[2], str(natoms))
                    if 'Max_Number_Events' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.max_events)
                    if 'Type_of_Events' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.type_events)
                    if 'Radius_Initial_Deformation' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.radius)
                    if 'Central_Atom' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.centre_atm)
                    if 'Eigenvalue_Threshold' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.eigen_thresh)
                    if 'Exit_Force_Threshold' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.force_thresh)
                    if 'Increment_Size' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.inc_size)
                    if 'Initial_Step_Size' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.init_step)
                    if 'sym_break_dist' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.sym_break)
                    if 'Max_Perp_Moves_Basin' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.max_perp_basin)
                    if 'Min_Number_KSteps' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.min_k)
                    if 'Basin_Factor' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.basin_factor)
                    if 'Lanczos_SCLoop' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.lanczos)
                    if 'Activation_MaxIter' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.activ_maxiter)
                    if 'delta_threshold' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.delta_thresh)
                    if 'Max_Perp_Moves_Activ' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.max_perp_activ)
                    if 'Force_Threshold_Perp_Rel' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.force_thresh_perp)
                    if 'Max_Iter_Basin' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.basin_maxiter)
                    if 'Write_xyz' in line:
                        line = line.replace(line.split('#')[0].split()[2], '.' + args.write_xyz + '.')
                    if 'Strategy_of_Search' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.search_strat)

                    m.write(line)

def create_refconfig():

    for each_file, directory in my_dict.iteritems():
        with open(each_file) as f: 
            with open (join(directory, refconfig), 'w+') as m:
                m.write(' run_id:         1000\n'  +  ' total_energy:   0\n')     
                for line in f:
                    if (re.findall(r'\s\d\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)):
                        coords = re.findall(r'\s\d\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)[0]
                        m.write(coords)

def create_filecounter():

    for directory in my_dict.itervalues():
        with open(join(directory,filecounter),'w+') as f:
            f.write("Counter:      1000")


def submitting_scripts():

    for each_file, directory in my_dict.iteritems():
        print("Creating submission files for " + directory)

        call(['sbatch', '-J', each_file.split('.')[0] + '_artgau', join(directory, submission_script)], shell = False) 



if __name__ == '__main__':

    create_directory()
    create_submission_file()
    create_exec_script()
    create_header_script()
    zipped_atomic_list = calculate_atoms()
    create_gauss_art(zipped_atomic_list)
    create_refconfig()
    create_filecounter()
    submitting_scripts()

