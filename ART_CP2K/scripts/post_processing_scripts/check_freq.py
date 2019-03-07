import numpy as np
import re
import argparse

parser = argparse.ArgumentParser(description = 'Add distance tolerance')
parser.add_argument('-f','--log_file',
        help = 'log file to test e.g. min1000.log')
parser.add_argument('-type','--type', required = True, 
        help = 'type of files to test')
args = parser.parse_args()

logfile_to_test = args.log_file

freq_list = []

def analyze_freq():

    with open (logfile_to_test) as f:
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
        init_coords = False
        new_line = ''
        for line in f:
            if "Input orientation" in line:
                init_coords = True
            if (init_coords):
                if (re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)):
                    new_line = new_line + re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)[0]
            if "Distance matrix" in line:
                break

    print(new_line)

    initial_coords = []
    for each_line in new_line.splitlines():
        init_coords_x = float(each_line.split()[0])
        init_coords_y = float(each_line.split()[1])
        init_coords_z = float(each_line.split()[2])
        
        initial_coords.append((init_coords_x, init_coords_y, init_coords_z))

    return initial_coords

def get_optimized_coords():

    with open (logfile_to_test) as f:
        for line_number, each_line in enumerate(f):
            if "Input orientation" in each_line:
                opt_coord_start_point = line_number
            if "Distance matrix" in each_line:
                opt_coord_end_point = line_number

    with open (logfile_to_test) as m:
        opt_coords = False
        new_line = ''
        for number, line in enumerate(m):
            if "Input orientation" in line and number == opt_coord_start_point:
                opt_coords = True
            if (opt_coords):
                if (re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)):
                    new_line = new_line + re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)[0]
            if "Distance matrix" in line and number == opt_coord_end_point:
                break

    print(new_line)

    optimized_coords = []
    for each_line in new_line.splitlines():
        opt_coords_x = float(each_line.split()[0])
        opt_coords_y = float(each_line.split()[1])
        opt_coords_z = float(each_line.split()[2])
        
        optimized_coords.append((opt_coords_x, opt_coords_y, opt_coords_z))

    return optimized_coords

def compare_coords(initial_coords, optimized_coords):
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

    if (np.allclose(dist_matrix_map['ART'], dist_matrix_map['Gaussian'], atol = 0.1)):
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


