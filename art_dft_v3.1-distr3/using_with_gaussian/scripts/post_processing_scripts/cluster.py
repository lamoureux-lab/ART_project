#TODO This program does this...

import glob
import numpy as np
import os.path


class Clusters:
    """
    Used to assign min_file results to different cluster groups
    """
    def __init__(self, filetype):
        #should be 'min' or 'sad'
        if(filetype is 'min' or filetype is 'sad'):
            self.filetype = filetype
        else:
            print 'Filetype set must be min or sad'
            exit()

        self.art_files = self.get_all_art_result_filenames()
        self.map_to_cluster = {}
        self.set_output_file('clusters.log')
        self.cluster = {}

    def set_output_file(self, output_file_name):
        self.output_file = output_file_name

    def make_json_list(self):
        """
        Creates
        :param map_to_cluster:
        :return:
        """
        json_clusters = {}
        for key, value in self.map_to_cluster.iteritems():
            if not value in json_clusters:
                json_clusters[value] = []

            json_clusters[value].append(key)

        print json_clusters


    def get_all_art_result_filenames(self):
        """
        This function will search the current working directory for all file names
        of the set filetype
        :return:
        """
        #Excludes .xyz files from cluster parsing
        return [fn for fn in glob.glob(self.filetype + '1' + '*')
                 if not os.path.splitext(fn)[1] == '.xyz']

    def calculate_cluster_map(self, file_list):
        dict = {}  # Initialize a Dictionary
        for file_name in file_list:
            coords = []
            with open(file_name) as f:
                for i in range(0, 3):
                    _ = f.readline()  # Skipping top 3 lines before the coordinates
                for line in f:
                    coord = line.split()
                    x = float(coord[1])
                    y = float(coord[2])
                    z = float(coord[3])
                    coords.append((x, y, z))
            dict[file_name] = coords  # Added coordinates to the dictionary

        for key in sorted(dict.keys()):
            mat = np.zeros((len(dict[key]), len(dict[key])))
            for n, (x0, y0, z0) in enumerate(dict[key]):
                for m, (x1, y1, z1) in enumerate(dict[key]):
                    dist = np.linalg.norm(np.array([x0, y0, z0]) - np.array([x1, y1, z1]))
                    mat[n, m] = dist
            different_from_all = True
            k_to_map = key

            for k in sorted(self.cluster.keys()):
                if np.allclose(mat, self.cluster[k], atol=1e-2):
                    different_from_all = False
                    k_to_map = k
                    break
            if different_from_all:
                self.cluster[key] = mat
            self.map_to_cluster[key] = k_to_map

        # return map_to_cluster

clusters = Clusters('min')
clusters.get_all_art_result_filenames()
clusters.calculate_cluster_map(clusters.art_files)
clusters.make_json_list()