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

reading_coords = False
reading_forces = False

with open(gaussian_log_file) as log:
    with open(gaussian_output_file,'w+') as g:
        for line in log:
            if "Input orientation" in line:
                reading_coords = True
                g.write("outcoor:\n")
            if "Distance matrix" in line:
                reading_coords = False
            if reading_coords:
                if (re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)):
                    coords = re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)[0]
                    g.write(coords)

            if "E(" in line:
               g.write("energy:\n")
               g.write(str(float(line.split()[4])*27.2113838668) + '\n')
               g.write("forces:\n")

            if "Forces (Hartrees/Bohr)" in line:
                reading_forces = True
            if "Cartesian Forces" in line:
                reading_forces = False

            if reading_forces: 
                if (re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)):
                    forces = re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)[0]
                    g.write(forces)
                
            if 'Stationary' in line:
                g.write('Stationary:')
                break






