import argparse
import cluster


parser = argparse.ArgumentParser(description = 'check member files of each cluster')
parser.add_argument('-m','--min_files', nargs = '*', help = 'files_to_test')
parser.add_argument('-tol','--dist_tol', default = 0.1, help = 'distance_tolerance')
args = parser.parse_args()

files_to_test = args.min_files
tolerance = float(args.dist_tol)

mapping = cluster.calculate_cluster_map(files_to_test, tolerance)
file_dict = cluster.order_cluster_mapping(mapping)


rep_name = ''    #name of representative file (representative of the cluster)
member_name = '' #name of member file
for key, value in file_dict.items():
    rep_name = key + '_results.txt'
    for i in range(len(value)):
        if not value[i] == key:
            member_name = value[i] + '_results.txt' 
            with open(rep_name, 'r') as rep:
                for line in rep:
                    if ('Success!' in line) or ('Error!' in line):
                        with open(member_name, 'w+') as m:
                            m.write(line + '\n' + 'Belongs to ' + key) 
