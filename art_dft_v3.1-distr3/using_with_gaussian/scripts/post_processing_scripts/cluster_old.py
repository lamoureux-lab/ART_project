#TODO This program does this...

import glob
import numpy as np

import os.path
from os import makedirs
from os.path import isdir
from shutil import copy


def get_art_files(filetype):
    #Excludes .xyz files from cluster parsing
    return [fn for fn in glob.glob(filetype + '*')
             if not os.path.splitext(fn)[1] == '.xyz' and not os.path.splitext(fn)[1] == '.log']

def make_json_list(map):

        json_clusters = {}
        for key, value in map.iteritems():
            if not value in json_clusters:
                json_clusters[value] = []

            json_clusters[value].append(key)

        return json_clusters


def calculate_cluster_map(file_list, tolerance):
    my_dict = {}  # Initialize a Dictionary
    for min_file in file_list:
        coords = []
        with open(min_file) as f:
            for i in range(0, 3):
                _ = f.readline()  # Skipping top 3 lines before the coordinates
            for line in f:
                coord = line.split()
                x = float(coord[1])
                y = float(coord[2])
                z = float(coord[3])
                coords.append((x, y, z))
        my_dict[min_file] = coords  # Added coordinates to the dictionary

    cluster = {}
    map_to_cluster = {}
    for key in sorted(my_dict.keys()):
        mat = np.zeros((len(my_dict[key]), len(my_dict[key])))
        for n, (x0, y0, z0) in enumerate(my_dict[key]):
            for m, (x1, y1, z1) in enumerate(my_dict[key]):
                dist = np.linalg.norm(np.array([x0, y0, z0]) - np.array([x1, y1, z1]))
                mat[n, m] = dist
        different_from_all = True
        k_to_map = key

        for k in sorted(cluster.keys()):
            if (np.allclose(mat, cluster[k], atol= tolerance)):
                different_from_all = False
                k_to_map = k
                break
        if (different_from_all):
            cluster[key] = mat
        #print(k_to_map)
        map_to_cluster[key] = k_to_map

    return map_to_cluster


filetype = 'min1'
tolerance = 0.01
files = get_art_files(filetype)
map_to_cluster = calculate_cluster_map(files, tolerance)
make_json_list(map_to_cluster)

