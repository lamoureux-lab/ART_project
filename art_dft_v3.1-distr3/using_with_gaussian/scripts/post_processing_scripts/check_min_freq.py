import sys
import re
import numpy as np

reading_init = False
reading_init_count = 0
init_text = ''
reading_opt = False
opt_text = ''
frequency = []
j = 0
for line in sys.stdin:
    if line.startswith(' 1\\1\\GINC-N069\\FOpt'):
        reading_opt = True
    elif line.startswith('                          Input orientation:'):
        reading_init = True
        reading_init_count += 1
    elif line.startswith('                    Distance matrix (angstroms):'):
        reading_init = False
    elif line.startswith(" Frequencies"):
    	j = j + 1
        frequency.append(line)
        print(line)

    if reading_opt:
        opt_text = opt_text + line
    if '\\@' in line:
        reading_opt = False

    pattern = re.compile('^\s+\d+\s+\d+\s+\d+\s+\S+\s+\S+\s+\S+$', flags=re.MULTILINE)
    if reading_init and (reading_init_count == 1):
        if pattern.match(line):
            init_text += line

s = init_text.split('\n')
init_coordinates = []
for line in s:
    if line == '':
        break
    coord = line.split()
    x_init = float(coord[3])
    y_init = float(coord[4])
    z_init = float(coord[5])
    init_coordinates.append((x_init,y_init,z_init))

opt_text = re.sub(r"^ ", "", opt_text, flags=re.MULTILINE)
opt_text = re.sub(r"\n", "", opt_text, flags=re.MULTILINE)

coord_vector = opt_text.split('\\\\')[3].split('\\')

opt_coordinates = []
for i in range(1,len(coord_vector)):
    coord_opt = coord_vector[i].split(',')
    x_opt = float(coord_opt[1])
    y_opt = float(coord_opt[2])
    z_opt = float(coord_opt[3])
    opt_coordinates.append((x_opt,y_opt,z_opt))

for line in frequency:
    freq = line.split()
    check = float(freq[2])
    if check < 0:
    	print("Optimization failed \n")
    	break

if j == 0:
	print("Failure: No frequency information found! \n")
elif j > 0:
	if check > 0:
		print ("Optimization Successful \n")

print('Initial Coordinates: \n')
print(init_coordinates) 
print('\n')

print('Optimized Coordinates: \n')
print(opt_coordinates) 
print('\n')

my_dict = {'initial' : init_coordinates, 'optimized' : opt_coordinates}

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
            if (np.allclose(mat, cluster[k], atol=1e-2)):
                different_from_all = False
                k_to_map = k
                break
        if(different_from_all):
            cluster[key] = mat
        map_to_cluster[key] = k_to_map

print("Clusters: \n")
print(map_to_cluster)
print('\n')

if(different_from_all):
    print('Error! ART and Gaussian coordinates are NOT similar \n')
else:
    print('Success! ART and Gaussian coordinates are similar \n')

 