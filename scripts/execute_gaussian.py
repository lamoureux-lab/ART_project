from subprocess import call
import re

gaussian_input_file = 'art2gaussian.inp'

gaussian_log_file = 'art2gaussian.log'

gaussian_output_file = 'gaussian2art'

with open(gaussian_input_file) as f:
    with open('combined_gaussian_input', 'a+') as c:
        for line in f:
            c.write(line)

call(['g09', gaussian_input_file], shell = False)


def reading_coords():
    stringed_coords = ''
    reading_coords = False
    with open(gaussian_log_file) as log:
        natoms = 0
        for line in log:
            if "Input orientation" in line:
                coords = []
                reading_coords = True
            if "Distance matrix" in line:
                reading_coords = False
            if reading_coords:
                if (re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)):
                    coords.append([float(line.split()[3]), float(line.split()[4]), float(line.split()[5])])
                    natoms += 1
            if 'Stationary' in line:
                break

        sum_x = 0
        sum_y = 0
        sum_z = 0
        for i in range(len(coords)):
            sum_x += (coords[i][0])
            sum_y += (coords[i][1])
            sum_z += (coords[i][2])

        centroid_x = (sum_x/natoms)
        centroid_y = (sum_y/natoms)
        centroid_z = (sum_z/natoms)

        for i in range(len(coords)):
            (coords[i][0]) -= centroid_x
            (coords[i][1]) -= centroid_y
            (coords[i][2]) -= centroid_z
        
        for each_list in coords:
            stringed_coords += (' ').join(str(x) for x in each_list) + '\n'
    return stringed_coords

def reading_energy():
    with open(gaussian_log_file) as log:
        for line in log:
            if "E(" in line:
                energy_string = ('energy:' + '\n' + str(float(line.split()[4])*27.2113838668) + '\n')
            if 'Stationary' in line:
                break
    return energy_string

def reading_forces():
    stringed_forces = ''
    reading_forces = False
    with open(gaussian_log_file) as log:
        for line in log:
            if "Forces (Hartrees/Bohr)" in line:
                forces = []
                reading_forces = True
            if "Cartesian Forces" in line:
                reading_forces = False

            if reading_forces: 
                if (re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)):
                    forces.append([float(line.split()[2]), float(line.split()[3]), float(line.split()[4])])

            if 'Stationary' in line:
                break
        
        for each_list in forces:
            stringed_forces += (' ').join(str(x) for x in each_list) + '\n'
    return stringed_forces               
        
if __name__ == "__main__":
    
    coords = reading_coords()    
    energy = reading_energy()    
    forces = reading_forces()    
    
    with open(gaussian_output_file, 'w+') as out:
        out.write('outcoor:' + '\n' + coords + energy + '\n' + 'forces:' + '\n' + forces + '\n' + 'Stationary:')


