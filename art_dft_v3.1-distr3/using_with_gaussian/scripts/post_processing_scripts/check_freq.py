import sys
frequency = []
for line in sys.stdin:
    if line.startswith(" Frequencies"):
        frequency.append(line)
        print(line)

for line in frequency:
    freq = line.split()
    check = float(freq[2])
    if check < 0:
    	break
    	print("Optimization failed")

if check > 0:
	print ("Optimization Successful")