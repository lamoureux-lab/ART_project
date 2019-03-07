import networkx as nx
import argparse
import numpy as np
import matplotlib.pyplot as plt
import re

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--files', help = 'add files for networking', nargs = '+')
args = parser.parse_args()

event_file = 'checked_events.list'

event_list = []
all_files = []
with open(event_file) as e:
    for line in e:
        init_min = line.split()[0]
        saddle = line.split()[1]
        final_min = line.split()[2]
        event_list.append((init_min, saddle, final_min))

for each_event in event_list:
    for each_file in each_event:
        all_files.append(each_file)

coords_dict = {}
centroid_x = {}
centroid_y = {}
centroid_z = {}
for each_file in all_files:
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
for each_file in all_files:
    member_list = []
    for every_other_file in all_files:
        cov = np.matmul(np.transpose(coords_dict[each_file]), coords_dict[every_other_file])
        u, s, v = (np.linalg.svd(cov))
        if(np.linalg.det(u)*np.linalg.det(v) < 0):
            v[:,-1]= -v[:,-1]
        rot = (np.matmul(u, v))
        
        rotated_coords = np.transpose(np.matmul(rot, np.transpose(coords_dict[every_other_file])))

        atoms_align_well = 0 
        for i in range(len(coords)):
            deviation_each_atom = np.sqrt(((rotated_coords[i][0] - coords_dict[each_file][i][0])**2    
            + (rotated_coords[i][1] - coords_dict[each_file][i][1])**2  
            + (rotated_coords[i][2] - coords_dict[each_file][i][2])**2))

            if deviation_each_atom < 0.1:
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

for (init_min, saddle, final_min) in event_list:
    for key, value in cluster_refined.items():
        if init_min in value:
            init_node = key
        if final_min in value:
            final_node = key
        if saddle in value:
            saddle_edge = key
    G.add_nodes_from([init_node, final_node])
    G.add_edge(init_node, final_node, title = saddle_edge)
    edge_name= nx.get_edge_attributes(G, 'title')

pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels = True)
nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_name)
plt.show()
