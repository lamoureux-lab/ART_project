import re
from os.path import join, dirname, relpath, isfile
from subprocess import call
import argparse
import cluster

submission_script_template = join(dirname(relpath(__file__)), '../gauss_graham.sub')
check_frequency            = join(dirname(relpath(__file__)), 'check_freq.py')
check_cluster_members      = join(dirname(relpath(__file__)), 'check_cluster_members.py')

parser = argparse.ArgumentParser(description = 'Create input and submission files and do gaussian optimization and then check the optimization ')

parser.add_argument('-f', '--filenames', nargs='+',
        help='specific input files to submit from project directory (e.g., min1000')
parser.add_argument('-i','--gaussian_input_template', 
        help='input file to extract the method, basis set, charge-multiplicity and coordinates from')
parser.add_argument('-j','--job_name', default = 'check',
        help = 'Name of the job')
parser.add_argument('-w','--wall_time', default = '1:00:00',
        help = 'Wall time in hh:mm:ss')
parser.add_argument('-mem','--memory', default = '2000MB',
        help = 'Memory in MB e.g. 2000MB')
parser.add_argument('-np','--nodes_proc', default = '4',
        help = 'number of nodes and processors, e.g. 4')
parser.add_argument('-tol','--tol_def_clusters', type = float, default = 0.1,
        help = 'Tolerance for defining clusters eg. 0.01')
parser.add_argument('-c','--check_one_per_cluster', action = 'store_true',
        help = 'Option to automatically check one file per cluster')

args = parser.parse_args()

file_type = re.split('(\d+)', args.filenames[0])[0]

input_file_template = args.gaussian_input_template
wt                  = args.wall_time
jn                  = args.job_name
np                  = args.nodes_proc
me                  = args.memory
tol_define_clusters = args.tol_def_clusters

def get_gaussian_header(input_file_template):
    count_link0_route = 0
    link0 = ''
    route = ''
    with open(input_file_template) as f:
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
                if file_type == 'min':
                    line = line.replace('force', 'opt=(maxcycles=400) freq ')
                elif file_type == 'sad':
                    line = line.replace('force', 'opt=(QST3, maxcycles=400) freq ')

                route = route + line

    return [link0, route, title, charge_mult]


def min_input_file(min_file, header_elements):
    with open (min_file) as f:
        coords = ''
        for line in f:
            if (re.findall(r'\S+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)):
                coords = coords + line

    min_input_file = min_file + '.inp'

    with open(min_input_file, 'w') as f:
        f.write(header_elements[0] + header_elements[1] + '\n' + min_file + '\n\n' + header_elements[3] + coords + '\n')

    return min_input_file

def sad_input_file(sad_file, header_elements):
    event = ()
    index = re.split('(\d+)', sad_file)[1]
    initial_min = 'min' + str(int(index) - 1)
    final_min = 'min' + index
    event = (initial_min, final_min, sad_file)
    
    sad_input_file = sad_file + '.inp'
    member_coords = {}
    for member_file in event:
        with open(member_file) as f:
            coords = ''
            for line in f:
                if (re.findall(r'\S+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)):
                    coords = coords + line
        member_coords[member_file] = coords

    with open(sad_input_file, 'w') as m:
        m.write(header_elements[0] + header_elements[1] + '\n' + 
                initial_min + '\n\n' + header_elements[3] + member_coords[initial_min] + '\n' + 
                final_min + '\n\n' + header_elements[3] + member_coords[final_min] + '\n' + 
                sad_file + '\n\n' + header_elements[3] + member_coords[sad_file] + '\n')

    return sad_input_file

def write_submission_script(filename, file_string):

    submission_script = filename + '.sub'

    call_gaussian      = 'g16' + ' ' + filename + '.inp' + '\n' 
    call_check_freq    = 'python ' + check_frequency + ' -f ' + filename + '.log' + ' -type ' + 'check_' + file_type + ' > ' + filename + '_results.txt' + '\n'  
    call_check_members = 'python ' + check_cluster_members + ' -f ' + file_string + '-tol ' + str(tol_define_clusters)

    with open(submission_script_template) as f:
        with open(submission_script,'w') as m:
            new_line = ''
            for line in f: 
                new_line = new_line + line 

            if args.check_one_per_cluster:
                m.write(new_line + call_gaussian + call_check_freq + call_check_members)

            else:
                m.write(new_line + call_gaussian + call_check_freq)

    return submission_script


def mapping_cluster_files(): 
     cluster_map = cluster.mapping(files_to_test, tol_define_clusters)
     return cluster_map


def string_of_files():
    file_string = ''
    for each_file in args.filenames:
        file_string = file_string + each_file + ' '
    return file_string


if __name__ == '__main__':

    header_elements = get_gaussian_header(input_file_template)
    files_to_test = args.filenames
    file_string = string_of_files()

    if args.check_one_per_cluster:
        cluster_map = mapping_cluster_files()
        files_to_test = []
        for key in cluster_map.keys():
            files_to_test.append(key)
    
    for filename in files_to_test:
        if not isfile(filename + '.log'):
            submission_script = write_submission_script(filename, file_string)
            if file_type == 'min':
                min_input_file(filename, header_elements)
            if file_type == 'sad':
                sad_input_file(filename, header_elements)
            print("Running jobs for ", filename)
        #    call(['sbatch','-J ' + jn + '_' + filename, '-t' + wt,'-c' + np, submission_script], shell=False)
      


