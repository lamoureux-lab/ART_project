import os
from os.path import join, dirname, relpath
from subprocess import call
import argparse
import random
import sys
import cluster


submission_script_template = join(dirname(relpath(__file__)), '../submission_script.sub')
checking_frequency_script =  join(dirname(relpath(__file__)), 'check_min_freq.py')
checking_cluster_members =  join(dirname(relpath(__file__)), 'check_cluster_members.py')



parser = argparse.ArgumentParser(description = 'Create input and submission files and do gaussian optimization and then check the optimization ')
parser.add_argument('-m', '--min_files', nargs='*',
        help='specific input files to submit from project directory (e.g., min1000')
parser.add_argument('-i','--gaussian_input_template', help='input file to extract the method, basis set, charge-multiplicity and coordinates from')
parser.add_argument('-j','--job_name', default = 'check',
        help = 'Name of the job')
parser.add_argument('-w','--wall_time', default = '1:00:00',
        help = 'Wall time in hh:mm:ss')
parser.add_argument('-mem','--memory', default = '2000MB',
        help = 'Memory in MB e.g. 2000MB')
parser.add_argument('-np','--nodes_proc', default = '4',
        help = 'number of nodes and processors, e.g. 4')
parser.add_argument('-tol1','--dist_tol', type = float, default = 0.1,
        help = 'Distance tolerance value eg. 0.01')
parser.add_argument('-tol2','--dist_tol2', type = float, default = 0.1,
        help = 'Distance tolerance value eg. 0.01')
parser.add_argument('-c','--check_one_per_cluster', action = 'store_true',
        help = 'Option to automatically check one file per cluster')


args = parser.parse_args()


input_file_template = args.gaussian_input_template
wt = args.wall_time
jn = args.job_name
np = args.nodes_proc
me = args.memory
tolerance = args.dist_tol
dist_tolerance = args.dist_tol2


def get_link0_route_title_charge_mult(input_file_template):
    count_link0_route = 0
    with open(input_file_template) as f:
        link0_route_title_charge_mult = ''
        for line_number, line in enumerate(f,1):
            if line.startswith('%') or line.startswith('#'):
                count_link0_route = count_link0_route + 1
            if line.startswith('#'):
                if 'opt' in line and 'freq' not in line:
                    line = line.replace('opt', 'opt=(maxcycles=400) freq ')
                elif 'opt' in line and 'freq' in line:
                    line = line.replace('opt', ' opt=(maxcycles=400) ')
                elif 'opt' not in line and 'freq' in line:
                    line = line.replace('freq', 'opt=(maxcycles=400) freq ')
                elif 'opt' not in line and 'freq' not in line:
                    line = line.replace('#', '# opt=(maxcycles=400) freq ')
            link0_route_title_charge_mult = link0_route_title_charge_mult + line
            if line_number == count_link0_route + 4: #After link0_route there is 'blankline', 'title', 'blankline', 'charge_multiplicity'
                break
        return link0_route_title_charge_mult


def get_atomic_coordinates(min_file):
    with open (min_file) as f:
        atomic_coords = ''
        for line_number, line in enumerate(f):
            if line_number > 2:
                atomic_coords = atomic_coords + line
        return atomic_coords


def create_gaussian_input_file(min_file):
    gaussian_input_file = min_file + '.inp'
    with open(gaussian_input_file, 'w') as f:
        f.write(link0_route_title_charge_mult + coords + '\n\n')

    return gaussian_input_file


def read_submission_template(submission_script_tempalte):
    with open(submission_script_tempalte) as f:
        new_line = ''
        for line in f:
            new_line = new_line + line
    return new_line


def write_submission_script(filename, submission_script_contents, dist_tolerance, file_string, check_per_cluster):
    new_line = ''
    submission_script = filename + '.sub'
    with open(submission_script,'w') as f:
        for line in submission_script_contents.splitlines():
            if line.startswith('g16'):
                line = line.replace('g16', 'g16 ' + filename + '.inp')
            elif line.startswith('python'):
                if check_per_cluster:
                    line = line.replace('python', 'python '+checking_frequency_script+' -tol2 '+ str(dist_tolerance)+' < '+min_file+'.log'+' > '+min_file+'_results.txt' + 
                                        '\n' + 'python '+ checking_cluster_members + ' -m ' + file_string + '-tol ' + str(args.dist_tol))
                if not check_per_cluster:
                    line = line.replace('python', 'python '+checking_frequency_script+' -tol2 '+ str(dist_tolerance)+' < '+min_file+'.log'+' > '+min_file+'_results.txt')
            
            new_line = new_line + line + '\n'
        f.write(new_line)
        return submission_script


def mapping_cluster_files(): 
     mapping = cluster.calculate_cluster_map(files_to_test, tolerance)
     file_dict = cluster.order_cluster_mapping(mapping)
     return(file_dict)


def string_of_files():
    file_string = ''
    for each_file in args.min_files:
        file_string = file_string + each_file + ' '
    return file_string


if __name__ == '__main__':

    link0_route_title_charge_mult = get_link0_route_title_charge_mult(input_file_template)
    submission_script_contents = read_submission_template(submission_script_template)



    files_to_test = args.min_files

    file_string = string_of_files()

    check_per_cluster = False

    if args.check_one_per_cluster:

        check_per_cluster = True

        file_dict = mapping_cluster_files()
        
        files_to_test = []
        for key in file_dict.keys():
            files_to_test.append(key)
    
  
    for min_file in files_to_test:
        if not os.path.isfile(min_file + '.inp'):
            submission_script = write_submission_script(min_file, submission_script_contents, dist_tolerance, file_string, check_per_cluster)
            coords = get_atomic_coordinates(min_file)
            create_gaussian_input_file(min_file)
            call(['sbatch','-J ' + jn + '_' + min_file, '-t' + wt,'-c' + np, submission_script], shell=False)
      


