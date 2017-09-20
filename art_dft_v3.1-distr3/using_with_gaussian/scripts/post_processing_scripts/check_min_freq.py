import sys
import re

reading_init = False
reading_init_count = 0
init_text = ''
reading_opt = False
opt_text = ''
frequency = []
j = 0
for line in sys.stdin:
    if line.startswith(' 1\\1\\GINC-N069\\FOpt'):
        reading_opt = True
    elif line.startswith('                          Input orientation:'):
        reading_init = True
        reading_init_count += 1
    elif line.startswith('                    Distance matrix (angstroms):'):
        reading_init = False
    elif line.startswith(" Frequencies"):
    	j = j + 1
        frequency.append(line)
        print(line)

    if reading_opt:
        opt_text = opt_text + line
    if '\\@' in line:
        reading_opt = False

#   pattern = re.compile('^\s+\d+\s+\d+\s+\d+\s+\w+\s+\w+\s+\w+$', flags=re.MULTILINE)
    pattern = re.compile('^\s+\d+\s+\d+\s+\d+\s+\S+\s+\S+\s+\S+$', flags=re.MULTILINE)
    if reading_init and (reading_init_count == 1):
        if pattern.match(line):
            init_text += line

print('INIT_TEXT')
print(init_text)
### TODO: Transform init_text into coords_init (same format as coords)

opt_text = re.sub(r"^ ", "", opt_text, flags=re.MULTILINE)
opt_text = re.sub(r"\n", "", opt_text, flags=re.MULTILINE)
coord_vector = opt_text.split('\\\\')[3].split('\\')
coords = []
for i in range(1,len(coord_vector)):
    coords.append(coord_vector[i].split(',')[1:])
print(coords)

### Compare coords and coords_init (with distance matrices)
### Return error if they are too different

for line in frequency:
    freq = line.split()
    check = float(freq[2])
    if check < 0:
    	print("Optimization failed")
    	break

if j == 0:
	print("Failure: No frequency information found!")
elif j > 0:
	if check > 0:
		print ("Optimization Successful")