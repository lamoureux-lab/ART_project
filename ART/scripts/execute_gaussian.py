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
    reading_coords = False
    with open(gaussian_log_file) as log:
        for line in log:
            if "Input orientation" in line:
                coords = ''
                reading_coords = True
            if "Distance matrix" in line:
                reading_coords = False
            if reading_coords:
                if (re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)):
                    coords += (re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line))[0]
            if 'Stationary' in line:
                break

    return coords

def reading_energy():
    with open(gaussian_log_file) as log:
        for line in log:
            if "E(" in line:
                energy = str(float(line.split()[4])*27.2113838668)
            if 'Stationary' in line:
                break
    return energy

def reading_forces():
    reading_forces = False
    with open(gaussian_log_file) as log:
        for line in log:
            if "Forces (Hartrees/Bohr)" in line:
                forces = ''
                reading_forces = True
            if "Cartesian Forces" in line:
                reading_forces = False
            if reading_forces: 
                if (re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)):
                    forces += (re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line))[0]
            if 'Stationary' in line:
                break
        
    return forces               
        
if __name__ == "__main__":
    
    coords = reading_coords()    
    energy = reading_energy()    
    forces = reading_forces()    
    
    with open(gaussian_output_file, 'w+') as out:
        out.write('outcoor:' + '\n' + coords + '\n' + 'energy:' + '\n' + energy + '\n\n' + 'forces:' + '\n' + forces) 


