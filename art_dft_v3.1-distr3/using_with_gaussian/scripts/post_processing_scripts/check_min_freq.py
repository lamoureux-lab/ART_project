import sys
frequency = []
j = 0
for line in sys.stdin:
    if line.startswith(" Frequencies"):
    	j = j + 1
        frequency.append(line)
        print(line)

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