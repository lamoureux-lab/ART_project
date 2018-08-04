import argparse
import cluster
from os.path import isfile


parser = argparse.ArgumentParser(description = 'check member files of each cluster')
parser.add_argument('-f','--files', nargs = '*', help = 'files_to_test')
parser.add_argument('-tol','--tol_define_clusters', default = 0.1, help = 'distance_tolerance')
args = parser.parse_args()

files_to_test = args.files
tolerance = float(args.tol_define_clusters)

cluster_map = cluster.mapping(files_to_test, tolerance)

rep_name = ''    #name of representative file (representative of the cluster)
member_name = '' #name of member file
for key, value in cluster_map.items():
    rep_name = key + '_results.txt'
    for i in range(len(value)):
        if not value[i] == key:
            member_name = value[i] + '_results.txt' 
            with open(rep_name, 'r') as rep:
                for line in rep:
                    if ('verdict' in line) :
                        if not isfile(member_name):
                            with open(member_name, 'w+') as m:
                                m.write(line + '\n' + 'Belongs to ' + key) 
