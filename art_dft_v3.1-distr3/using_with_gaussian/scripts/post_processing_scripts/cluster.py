#Cluster the min and sad files based on their distance matrices.

import numpy as np


def order_cluster_mapping(cluster_map):
    
    ordered_cluster_map = {}
    for key, value in sorted(cluster_map.items()): # Changed iteritems to items for python3
        if not value in ordered_cluster_map:
            ordered_cluster_map[value] = []
        ordered_cluster_map[value].append(key)
    return ordered_cluster_map


def calculate_cluster_map(file_list, tolerance):
    my_dict = {}  # Initialize a Dictionary
    for min_file in file_list:
        coords = []
        with open(min_file) as f:
            for line_number, line in enumerate(f,1):
                if line_number > 3:
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
        map_to_cluster[key] = k_to_map
    return map_to_cluster



