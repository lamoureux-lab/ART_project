import networkx as nx
import argparse
import numpy as np
import matplotlib.pyplot as plt
import re

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--files', help = 'add files for networking', nargs = '+')
args = parser.parse_args()

coords_dict = {}
centroid_x = {}
centroid_y = {}
centroid_z = {}
for each_file in args.files:
    with open(each_file) as f:
        coords = []
        natoms = 0
        for line in f:
            if (re.findall(r'\S+\s+[-]?\d+[.]\d+\s+[-]?\d+[.]\d+\s+[-]?\d+[.]\d+\s+',line)): 
                coords.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
                natoms = natoms + 1
        coords_dict[each_file] = coords

        centroid_x[each_file] = np.mean(coords_dict[each_file], axis = 0)[0]
        centroid_y[each_file] = np.mean(coords_dict[each_file], axis = 0)[1]
        centroid_z[each_file] = np.mean(coords_dict[each_file], axis = 0)[2]
        
        for i in range(len(coords)):
            coords_dict[each_file][i][0] = coords_dict[each_file][i][0] - centroid_x[each_file] 
            coords_dict[each_file][i][1] = coords_dict[each_file][i][1] - centroid_y[each_file]
            coords_dict[each_file][i][2] = coords_dict[each_file][i][2] - centroid_z[each_file]
cluster = {}

for each_file in sorted(args.files):
    member_list = []
    for every_other_file in args.files:
        cov = np.matmul(np.transpose(coords_dict[each_file]), coords_dict[every_other_file])
        u, s ,v = (np.linalg.svd(cov))
        v[:,-1]= -v[:,-1]
        if(np.linalg.det(u)*np.linalg.det(v) < 0):
            v[:,-1]= -v[:,-1]
        rot = (np.matmul(u, v))
        
        rotated_coords = np.transpose(np.matmul(rot, np.transpose(coords_dict[every_other_file])))

        atoms_align_well = 0 
        for i in range(len(coords)):
            rmsd_each_atom = np.sqrt(((rotated_coords[i][0] - coords_dict[each_file][i][0])**2    
            + (rotated_coords[i][1] - coords_dict[each_file][i][1])**2  
            + (rotated_coords[i][2] - coords_dict[each_file][i][2])**2)/natoms)

            if rmsd_each_atom < 0.05:
                atoms_align_well += 1

        sum_rmsd = 0
        for i in range(len(coords)):
            sum_rmsd += ((rotated_coords[i][0] - coords_dict[each_file][i][0])**2    
            + (rotated_coords[i][1] - coords_dict[each_file][i][1])**2  
            + (rotated_coords[i][2] - coords_dict[each_file][i][2])**2)

        rmsd = np.sqrt((sum_rmsd/natoms))
        if rmsd < 0.1 and atoms_align_well == natoms:
            member_list.append(every_other_file)
            cluster[each_file] = member_list

cluster_refined = {}

for rep, members in sorted(cluster.items()):
    if members not in cluster_refined.values():
        cluster_refined[rep] = members


print(cluster_refined)

G = nx.Graph()

G.add_nodes_from(cluster_refined.keys())

nx.draw(G, with_labels = True)

plt.show()
