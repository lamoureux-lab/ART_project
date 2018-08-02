import argparse
import re
import numpy as np

def mapping(files_to_test, tolerance):
    file_coords = {}
    for each_file in files_to_test:
        coords = []
        with open(each_file) as f:
            for line in f:
                if (re.findall(r'\s\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)):
                    coords.append((float(line.split()[1]), float(line.split()[2]), float(line.split()[3])))
            file_coords[each_file] = coords    

    file_matrix = {}
    for every_file, coordinates in file_coords.items():
        distance_matrix = []
        for each_atom in coordinates:
            for every_other_atom in coordinates:
                distance = ((each_atom[0] - every_other_atom[0])**2 + (each_atom[1] - every_other_atom[1])**2 + (each_atom[2] - every_other_atom[2])**2)**0.5
                distance_matrix.append(distance)
    
        file_matrix[every_file] = distance_matrix

    cluster = {}
    count = 0
    for every_file in sorted(file_matrix.keys()):
        member_list = []
        for every_other_file in sorted(file_matrix.keys()):
            if (np.allclose(file_matrix[every_file], file_matrix[every_other_file], atol = tolerance)):
                    member_list.append(every_other_file)
                    cluster[every_file] = member_list

        for member in member_list:
            for each_test_file in files_to_test:
                if member == each_test_file:
                    count = count + 1

        if count == len(files_to_test):
            break

    return cluster

