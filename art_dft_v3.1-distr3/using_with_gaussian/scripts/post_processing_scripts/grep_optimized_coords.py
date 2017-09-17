import itertools
import sys

key_phrase_1 = '    -- Stationary point found.'
key_phrase_2 = ' Center'
key_phrase_3 = '                    Distance matrix (angstroms):'

#Finds the starting positions of all the coordinates whether optimized or not and stores them in a list
def find_coordinate_position(input_file, keyphrase):
	phrase_list = []
	i = 0
	for line in input_file:
		i = i + 1
		if line.startswith(keyphrase):
			phrase_list.append(i)
	return phrase_list

#Finds the stationary point location in the file
def find_stationary_point(input_file, keyphrase):
	i = 0
	for line in input_file:
		i = i + 1
		if line.startswith(keyphrase):
			break
	return i

#Returns the start and end points of the optimized coords
def get_optimized_coords(coord_list, point):
	for element in coord_list:
		if element > point:
			break
	return element

#Slices the file at the start and end point of the optimized coords
def slice_coords(input_file, start, end):
	s = ''
	for line in itertools.islice(input_file, start_optimized_coord, end_optimized_coord):
		s = s + line
	return s

#Writes the sliced content (i.e. the optimized coordinates) into a new file
def write_coords(output_file):
	file_out.write(coords)




stationary_point = find_stationary_point(sys.stdin, key_phrase_1)
coord_start_list = find_coordinate_position(sys.stdin, key_phrase_2)
coord_end_list = find_coordinate_position(sys.stdin, key_phrase_3)
start_optimized_coord = get_optimized_coords(coord_start_list, stationary_point)
end_optimized_coord = get_optimized_coords(coord_end_list, start_optimized_coord)
slice_coords(sys.stdin, start_optimized_coord, end_optimized_coord)
write_coords(sys.stdout)

