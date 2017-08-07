import sys
frequency = []
for line in sys.stdin:
    if line.startswith(" Frequencies"):
        frequency.append(line)
        print(line)
i = 0
for line in frequency:
    freq = line.split()
    check = float(freq[2])
    if check < 0:
    	i = i + 1

if i == 1:
	print ("Optimization Successful")
else:
	print("Optimization failed")