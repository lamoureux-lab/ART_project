import argparse
import re
parser = argparse.ArgumentParser(description = 'extract coords')
parser.add_argument('-s', '--files', nargs = '+', help='enter sad files (e.g., sad1000)')
args = parser.parse_args()

input_files = args.files

events = []
for sad_file in sorted(input_files):
    index = re.split('(\d+)', sad_file)[1]
    initial_min = 'min' + str(int(index) - 1)
    events.append((initial_min, sad_file))

final_min = 'min' + index  # The very last min (because input_files are sorted, so index is maximum outside the loop)

with open(final_min) as f: 
    natoms = 0
    final_min_coords = ''
    for line in f:
        if (re.findall(r'\s\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)):
            natoms = natoms + 1
            final_min_coords = final_min_coords + line

new_line = ''
for each_event in events:
    for member_file in each_event:
        new_line = new_line + str(natoms) + '\n' + member_file + '\n' 
        with open(member_file) as m:
            for line in m:
                if (re.findall(r'\s\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+[-]?\d[.]\d+\s+',line)):
                    new_line = new_line + line

print(new_line)
print(str(natoms) + '\n' + final_min + '\n' + final_min_coords)

