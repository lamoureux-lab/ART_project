from numpy import matrix, linalg as LA 

parser = argparse.ArgumentParser(description = 'Calculate displacement vector')
parser.add_argument('-s', '--sad_art_files', nargs='*',
                    help='specific ART sad output files to submit from project directory (e.g., sad1001')

def get_dispacement_vector(sad_file, initial_min_file):
    with open (sad_file) as f:
        atomic_coords = ''
        for i in range (0,3):
            _ = f.readline()
        for line in f:
        	atomic_coords = atomic_coords + line
            s = line.split()
        	x = float(s[1])
        	y = float(s[2])
        	z = float(s[3])
        	m1 = matrix([[x],[y],[z]])

    with open (initial_min_file) as g:
        atomic_coords = ''
        for i in range (0,3):
            _ = g.readline()
        for line in f:
            atomic_coords = atomic_coords + line
        	s = line.split()
        	x = float(s[1])
        	y = float(s[2])
        	z = float(s[3])
        	m2 = matrix([[x],[y],[z]])

    dv = m2 - m1
    mag = LA.norm(v)
    dv_norm = dv/mag
    return dv_norm

sad_files = args.sad_art_files

if __name__ == '__main__':

	for sad_file in sad_files:
		displacement_vector = get_dispacement_vector(sad_file, 'min1000')

