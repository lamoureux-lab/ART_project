from os.path import dirname, join, relpath, isfile, exists
from os import makedirs, getcwd, chdir
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
parser.add_argument('-am', '--activ_max', help = 'max iterations during activation, e.g. 400', default ='400')
parser.add_argument('-d', '--delta_thresh', help = 'energy threshold during Lanczos e.g. 4.0', default ='4.0')
parser.add_argument('-pa', '--max_perp_activ', help = 'max perp steps during activation e.g 5', default ='5')
parser.add_argument('-fp', '--force_thresh_perp', help = 'force thresh perp relax, e.g. 0.5', default ='0.5')
parser.add_argument('-bm', '--basin_max', help = 'max iterations in basin, e.g. 20', default ='20')
parser.add_argument('-xyz', '--write_xyz', help = 'write or not xyz files -- true or false', default ='true')
parser.add_argument('-search', '--search_strat', help = 'search strategy, e.g. 0, 1 or 2', default ='0')
parser.add_argument('-read', '--read_from', help = 'read from a specific file, e.g. vector20', default ='xxx')

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

shared_log =                 '../shared_log/'

files_to_run =                args.gaussian_input

my_dict =                     {}

for each_file in files_to_run:
    my_dict[each_file] = re.split(r'[.]', each_file)[0]

def create_directory():
    for directory in my_dict.values():
        try:
            if not exists(directory):
                makedirs(directory)
        except OSError:
            print ('Error: Creating directory. ' +  directory)

def create_submission_file():
    for directory in my_dict.values():        
        with open(submission_script_template) as f:
            with open(join(directory, submission_script), 'w+') as m:
                new_line = ''
                for line in f:
                    if 'mem' in line:
                        line = line.replace(line.split('=')[1], args.mem + '\n')
                    if 'nodes' in line:
                        line = line.replace(line.split('=')[1], args.nodes + '\n')
                    if 'time' in line:
                        line = line.replace(line.split('=')[1], args.walltime + '\n')
                    if 'gres' in line:
                        line = line.replace(line.split('=')[1], 'gpu:' + args.gpus + '\n')
                    if 'cpus' in line:
                        line = line.replace(line.split('=')[1], args.cpus + '\n')
                    new_line = new_line + line
                calling_art = 'csh ' + gauss_art + ' > ' + 'output.log'
                m.write(new_line + calling_art)
               
def create_exec_script():
    for directory in my_dict.values():
        with open(exec_script_template) as f:
            with open(join(directory, exec_script), 'w+') as m:
                for line in f:
                    m.write(line)

def create_header_script():
    for each_file, directory in my_dict.items():
        with open(create_header_template) as f:
            with open(join(directory, create_header),'w+') as m:
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
                if (re.findall(r'\S+\s+[-]?\d+[.]\d+\s+[-]?\d+[.]\d+\s+[-]?\d+[.]\d+\s+',line)):
                    natoms = natoms + 1
            
            atom_list.append(natoms)

    zipped_atomic_list = zip(files_to_run, atom_list)
    return zipped_atomic_list

def read_vec_log():
    vec_log = 'shared_log' + '/' + args.read_from + '.log'
    if isfile(vec_log):
        min_count = 0
        sad_count = 0
        disp_count = 0
        natoms_read = 0
        with open(vec_log) as f:
            for row in f:
                if 'min' in row:
                    min_count = min_count + 1
                if 'sad' in row:
                    sad_count = sad_count + 1
                if 'Displacement' in row:
                    disp_count = disp_count + 1
                if (re.findall(r'\S+\s+[-]?\d+[.]\d+\s+[-]?\d+[.]\d+\s+[-]?\d+[.]\d+\s+',row)):
                    natoms_read = natoms_read + 1

        natoms_read = int(natoms_read / (min_count + sad_count + disp_count)) 

        min_sad_natoms_read = (min_count, sad_count, natoms_read)

    if not isfile(vec_log):
        min_sad_natoms_read = (0,0,0)

    return min_sad_natoms_read

def create_gauss_art(zipped_atomic_list, min_sad_natoms_read):
    for each_file, natoms in zipped_atomic_list:
        with open(gauss_art_template) as f:
            with open(join(each_file.split('.')[0], gauss_art), 'w+') as m:
                for line in f:                     
                    if 'NATOMS' in line:
                        line = line.replace(line.split('#')[0].split()[2], str(natoms), 1)
                    if 'Max_Number_Events' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.max_events, 1)
                    if 'Type_of_Events' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.type_events, 1)
                    if 'Radius_Initial_Deformation' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.radius, 1)
                    if 'Central_Atom' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.centre_atm, 1)
                    if 'Eigenvalue_Threshold' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.eigen_thresh, 1)
                    if 'Exit_Force_Threshold' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.force_thresh, 1)
                    if 'Increment_Size' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.inc_size, 1)
                    if 'Initial_Step_Size' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.init_step, 1)
                    if 'sym_break_dist' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.sym_break, 1)
                    if 'Max_Perp_Moves_Basin' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.max_perp_basin, 1)
                    if 'Min_Number_KSteps' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.min_k, 1)
                    if 'Basin_Factor' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.basin_factor, 1)
                    if 'Lanczos_SCLoop' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.lanczos, 1)
                    if 'Activation_MaxIter' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.activ_max, 1)
                    if 'delta_threshold' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.delta_thresh, 1)
                    if 'Max_Perp_Moves_Activ' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.max_perp_activ, 1)
                    if 'Force_Threshold_Perp_Rel' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.force_thresh_perp, 1)
                    if 'Max_Iter_Basin' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.basin_max, 1)
                    if 'Write_xyz' in line:
                        line = line.replace(line.split('#')[0].split()[2], '.' + args.write_xyz + '.', 1)
                    if 'Strategy_of_Search' in line:
                        line = line.replace(line.split('#')[0].split()[2], args.search_strat, 1)
                    if 'Shared_History_Filename' in line:
                        line = line.replace(line.split('#')[0].split()[2], shared_log + each_file.split('.')[0] + '.log', 1)
                    if 'Shared_History_To_Read' in line:
                        line = line.replace(line.split('#')[0].split()[2], shared_log + args.read_from + '.log', 1)
                    if 'nmin_read' in line:
                        line = line.replace(line.split('#')[0].split()[2], str(min_sad_natoms_read[0]), 1)
                    if 'nsad_read' in line:
                        line = line.replace(line.split('#')[0].split()[2], str(min_sad_natoms_read[1]), 1)
                    if 'natoms_read' in line:
                        line = line.replace(line.split('#')[0].split()[2], str(min_sad_natoms_read[2]), 1)

                    m.write(line)

def create_refconfig():
    for each_file, directory in my_dict.items():
        with open(each_file) as f: 
            with open (join(directory, refconfig), 'w+') as m:
                m.write(' run_id:         1000\n')     
                m.write(' total_energy:   0 \n')     
                for line in f:
                    if (re.findall(r'\S+\s+[-]?\d+[.]\d+\s+[-]?\d+[.]\d+\s+[-]?\d+[.]\d+\s+',line)): 
                        coords = re.findall(r'\S+\s+[-]?\d+[.]\d+\s+[-]?\d+[.]\d+\s+[-]?\d+[.]\d+\s+',line)[0]
                        m.write(coords)

def create_filecounter():
    for directory in my_dict.values():
        with open(join(directory, filecounter),'w+') as f:
            f.write("Counter:      1000")


def submitting_scripts():
    for each_file, directory in my_dict.items():
        print("Creating submission files for " + directory)
        work_directory = getcwd()
        chdir(join(work_directory, directory))
#        call(['sbatch', '-J', each_file.split('.')[0], submission_script], shell = False) 
        chdir(work_directory)


if __name__ == '__main__':

    create_directory()
    create_submission_file()
    create_exec_script()
    create_header_script()
    zipped_atomic_list = calculate_atoms()
    min_sad_natoms_read = read_vec_log()
    create_gauss_art(zipped_atomic_list, min_sad_natoms_read)
    create_refconfig()
    create_filecounter()
    submitting_scripts()

