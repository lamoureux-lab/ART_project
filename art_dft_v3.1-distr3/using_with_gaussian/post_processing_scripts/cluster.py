#TODO This program does this...

import glob
import numpy as np

import os.path

filetype = 'min1'
#Excludes .xyz files from cluster parsing
files = [fn for fn in glob.glob(filetype + '*')
         if not os.path.splitext(fn)[1] == '.xyz']

dict = {} #Initialize a Dictionary
for min in files:
    coords = []
    with open(min) as f:
            for i in range(0, 3):
                _ = f.readline() #Skipping top 3 lines before the coordinates
            for line in f:
                coord = line.split()
                x = float(coord[1])
                y = float(coord[2])
                z = float(coord[3])
                coords.append((x,y,z))
    dict[min] = coords #Added coordinates to the dictionary
cluster = {}
map_to_cluster = {}
for key in sorted(dict.keys()):
    mat = np.zeros( (len(dict[key]), len(dict[key])) )
    for n, (x0, y0, z0) in enumerate(dict[key]):
        for m, (x1, y1, z1) in enumerate(dict[key]):
            dist = np.linalg.norm(np.array([x0,y0,z0]) - np.array([x1,y1,z1]))
            mat[n,m] = dist
    different_from_all = True
    k_to_map = key
    for k in sorted(cluster.keys()):
        print(k)
        if (np.allclose(mat, cluster[k], atol=1e-2)):
            different_from_all = False
            k_to_map = k
            break
    if (different_from_all):
        cluster[key] = mat
    print(k_to_map)
    map_to_cluster[key] = k_to_map

print("Clusters:")
print(cluster.keys())
print(cluster)
print("Mappings:")
print(map_to_cluster)
