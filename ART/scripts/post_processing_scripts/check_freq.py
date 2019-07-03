import numpy as np
import re
import argparse

parser = argparse.ArgumentParser(description = 'Checks log files for frequencies and geometry optimizations.')
parser.add_argument('-f','--log_file',
        help = 'log file to test e.g. min1000.log')
parser.add_argument('-type','--type', required = True, 
        help = 'type of files to test, e.g. check_min for min files and check_sad for sad files')
parser.add_argument('-tol_opt','--tol_check_opt', type = float, default = 0.1, help = 'distance_tolerance to check optimization')
args = parser.parse_args()

logfile_to_test = args.log_file

freq_list = []

def analyze_freq():

    with open (logfile_to_test) as f:
        frequency_verdict = ''
        freq_count = 0
        for line in f:
            if "Frequencies" in line:
                freq_count += 1
                for i in range(2,len(line.split())):
                    freq_list.append(float(line.split()[i]))
                    print((line.split()[i]) + '\n')

    if freq_count == 0:
        frequency_verdict = "Failure: No frequency information found!\n"

    else:
        if args.type == 'check_min':
            for frequency in freq_list:
                if frequency < 0:
                    frequency_verdict = "Not a true minimum :(\n"
                    break
         
            if frequency > 0:
                frequency_verdict = "Indeed a true minimum :)\n"

        elif args.type == 'check_sad':
            negative_freq = 0
            for frequency in freq_list:
                if frequency < 0:
                    negative_freq += 1

            if negative_freq == 1:
                frequency_verdict = "Indeed a true first-order saddle :)\n"
            else:
                frequency_verdict = "Not a true saddle :( \n"

    print("Frequency verdict: ", frequency_verdict)


def get_initial_coords():

    with open (logfile_to_test) as f:
        init_coords_found = False
        init_coords = ''
        for line in f:
            if "Input orientation" in line:
                init_coords_found = True
            if (init_coords_found):
                if (re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)):
                    init_coords = init_coords + re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)[0]
            if "Distance matrix" in line:
                break

    print(init_coords)

    initial_coords = []
    for each_coord in init_coords.splitlines():
        init_coords_x = float(each_coord.split()[0])
        init_coords_y = float(each_coord.split()[1])
        init_coords_z = float(each_coord.split()[2])
        
        initial_coords.append((init_coords_x, init_coords_y, init_coords_z))

    return initial_coords

def get_optimized_coords():

    with open (logfile_to_test) as m:
        opt_coords_found = False
        opt_coords = ''
        for line in m:
            if "Input orientation" in line:
                opt_coords_found = True
                opt_coords = ''
            if "Distance matrix" in line:
                opt_coords_found = False
            if (opt_coords_found):
                if (re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)):
                    opt_coords = opt_coords + re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)[0]
    print(opt_coords)

    optimized_coords = []
    for each_coord in opt_coords.splitlines():
        opt_coords_x = float(each_coord.split()[0])
        opt_coords_y = float(each_coord.split()[1])
        opt_coords_z = float(each_coord.split()[2])
        
        optimized_coords.append((opt_coords_x, opt_coords_y, opt_coords_z))

    return optimized_coords

def compare_coords(initial_coords, optimized_coords):
    optimization_verdict = ''
    coords_map = {}
    coords_map['ART'] = initial_coords
    coords_map['Gaussian'] = optimized_coords
    
    dist_matrix_map = {}
    for program, coordinates in coords_map.items():
        distance_matrix = []
        for each_atom in coordinates:
            for every_other_atom in coordinates:
                distance = ((each_atom[0] - every_other_atom[0])**2 + (each_atom[1] - every_other_atom[1])**2 + (each_atom[2] - every_other_atom[2])**2)**0.5 
                distance_matrix.append(distance)

        dist_matrix_map[program] = distance_matrix

    if (np.allclose(dist_matrix_map['ART'], dist_matrix_map['Gaussian'], atol = args.tol_check_opt)):
        optimization_verdict = "ART structure is OK :)"
    else:
        optimization_verdict = "ART structure is not OK :("

    print("Optimization verdict: ", optimization_verdict)

if __name__ == '__main__':

    print("Optimization Results\n")
    print("Checking frequencies\n")
    analyze_freq()
    print("Checking structure\n")
    print("ART coordinates\n")
    initial_coords = get_initial_coords()
    print("Gaussian optimized coordinates\n")
    optimized_coords = get_optimized_coords()
    compare_coords(initial_coords, optimized_coords)


