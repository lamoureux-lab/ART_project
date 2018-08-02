import argparse
import os
from os.path import join, abspath, dirname, split
from os import listdir, makedirs, getcwd
from os.path import isfile, join, exists, splitext, relpath
import glob
from subprocess import call
import random
import sys
import re
import cluster

submission_script_template = join(dirname(relpath(__file__)), '../submission_script.sub')
checking_frequency_script =  join(dirname(relpath(__file__)), 'check_sad_freq.py')
checking_cluster_members =  join(dirname(relpath(__file__)), 'check_cluster_members.py')


parser = argparse.ArgumentParser(description = 'Create input and submission files and do gaussian optimization and then check the optimization ')
parser.add_argument('-s', '--sad_files', nargs='*',
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
    link0 = ''
    route = ''
    
    with open(input_file_template) as f:
        link0_route_title_charge_mult = ''
        for line_number, line in enumerate(f,1):
            
            if line.startswith('%') or line.startswith('#'):
                count_link0_route = count_link0_route + 1
            
            if line.startswith('%'):
                link0 = link0 + line  
            if line_number == count_link0_route + 2:
                title = line
            if line_number == count_link0_route + 4:
                charge_mult = line
           
            if line.startswith('#'):
                if 'opt' in line and 'freq' not in line:
                    line = line.replace('opt', 'opt=(QST3, maxcycles=400) freq ')
                elif 'opt' in line and 'freq' in line:
                    line = line.replace('opt', ' opt=(QST3, maxcycles=400) ')
                elif 'opt' not in line and 'freq' in line:
                    line = line.replace('freq', 'opt=(QST3, maxcycles=400) freq ')
                elif 'opt' not in line and 'freq' not in line:
                    line = line.replace('#', '# opt=(QST3, maxcycles=400) freq ')
                
                route = route + line

        gaussian_elements = [link0, route, title, charge_mult] 
        return gaussian_elements
             

def create_gaussian_input_file(sad_file):
    index = []
    events = []
    index.append(re.split('(\d+)', sad_file)[1])
    for each_index in index:
        prev_index = int(each_index) - 1
        initial_min = 'min' + str(prev_index)
        final_min = 'min' + (each_index)
    events.append((initial_min, final_min, sad_file))

    for each_event in events:
        new_line = '' 
        #print(link0_route_title_charge_mult)
        for member_file in each_event:
            with open(member_file) as f:
                for line_number, line in enumerate(f):
                    if line_number == 2:
                        line = '\n' + str(member_file) + '\n\n' + charge_mult
                        new_line = new_line + line
                    if line_number > 2:
                        new_line = new_line + line

        with open(sad_file+'.inp', 'w') as s:
            s.write(link0 + route + new_line + '\n')


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
                    line = line.replace('python', 'python '+checking_frequency_script+' -tol2 '+ str(dist_tolerance)+' < '+sad_file+'.log'+' > '+sad_file+'_results.txt' + 
                                        '\n' + 'python '+ checking_cluster_members + ' -m ' + file_string + '-tol ' + str(args.dist_tol))
                if not check_per_cluster:
                    line = line.replace('python', 'python '+checking_frequency_script+' -tol2 '+ str(dist_tolerance)+' < '+sad_file+'.log'+' > '+sad_file+'_results.txt')
            
            new_line = new_line + line + '\n'
        f.write(new_line)
        return submission_script


def mapping_cluster_files(): 
     mapping = cluster.calculate_cluster_map(files_to_test, tolerance)
     file_dict = cluster.order_cluster_mapping(mapping)
     return(file_dict)


def string_of_files():
    file_string = ''
    for each_file in args.sad_files:
        file_string = file_string + each_file + ' '
    return file_string


if __name__ == '__main__':
    
    link0 = get_link0_route_title_charge_mult(input_file_template)[0]
    route = get_link0_route_title_charge_mult(input_file_template)[1]
    title = get_link0_route_title_charge_mult(input_file_template)[2]
    charge_mult = get_link0_route_title_charge_mult(input_file_template)[3]
    
    submission_script_contents = read_submission_template(submission_script_template)



    files_to_test = args.sad_files

    file_string = string_of_files()

    check_per_cluster = False

    if args.check_one_per_cluster:

        check_per_cluster = True

        file_dict = mapping_cluster_files()
        
        files_to_test = []
        for key in file_dict.keys():
            files_to_test.append(key)
    
  
    for sad_file in files_to_test:
        if not os.path.isfile(sad_file + '.log'):
            submission_script = write_submission_script(sad_file, submission_script_contents, dist_tolerance, file_string, check_per_cluster)          
            create_gaussian_input_file(sad_file)
            call(['sbatch','-J ' + jn + '_' + sad_file, '-t' + wt,'-c' + np, submission_script], shell=False)
 


