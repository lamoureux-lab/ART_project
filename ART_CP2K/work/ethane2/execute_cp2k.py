from subprocess import call
import re

cp2k_log_file = 'cp2k_energy_forces'

call(['cp2k.popt', 'cp2k_input.inp'], shell = False)

def reading_energy():
    energy = ''
    with open (cp2k_log_file) as log:
        for line in log:
            if "ENERGY|" in line:
                energy = line.split(':')[1]
            elif "GEOMETRY OPTIMIZATION COMPLETED" in line:
                break
    return energy

def reading_coords():
    coords = ''
    with open (cp2k_log_file) as log:
        for line in log:
            if 'ATOMIC COORDINATES IN angstrom' in line:
                reading_coords = True
                coords = ''
            elif 'SCF PARAMETERS' in line:
                reading_coords = False
            if reading_coords:
                if (re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)):
                    coords += (re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line))[0]
     return coords

def reading_forces():
    reading_forces = False
    forces = ''
    with open(cp2k_log_file) as log:
        for line in log:
            if "ATOMIC FORCES in [a.u.]" in line:
                reading_forces = True
                forces = ''
            elif "SUM OF ATOMIC FORCES" in line:
                reading_forces = False
            if reading_forces: 
                if (re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)):
                    forces += (re.findall(r'[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line))[0]
            elif "GEOMETRY OPTIMIZATION COMPLETED" in line:
                break
    return forces               


if __name__ == "__main__":

    energy = reading_energy()

    forces = reading_forces()

    with open ("cp2k_output", 'w+') as f:
        f.write("coords: \n\n" + coords + "energy: \n\n" + energy + "forces: \n\n" + forces)


