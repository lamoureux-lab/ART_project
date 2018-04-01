class myDict(dict):       # Defined a class to append elements to the dictionary in the 'sorting' function, see "def sorting" below
    def __init__(self):
        self = dict()

    def add(self, key, value):
        self[key] = value 


def grep_activation_relaxation(preferred_event, non_preferred_event):
	with open("output.log") as f:
		occurence = False
		new_line = ''
		count = 0
		for line in f:
			count = count + 1
          		line = str(count) + '\t' + line
			if preferred_event in line:
				occurence = True
			if non_preferred_event in line:		
				occurence = False
			if occurence:
				new_line = new_line + line
                return new_line



def coords_before_K0():
    with open('output.log') as f:
        count = 0
        coords = False
        new_line = ''
        for line in f:
            count = count + 1
            line = str(count) + '\t' + line
            if "posa" in line:
                coords = True
            if "forces" in line:
                coords = False
            if coords:
                new_line = new_line + line
            if "K converged" in line:
                break
        return new_line



def remove_junk_lines(event_type, event_name):
	coords = False
	new_line = ''
	for line in event_type.splitlines():
		if "posa" in line:
			coords = True
			line = line.replace("posa", event_name)
		if "forces:" in line:
			coords = False
		if coords:
			line = line + '\n'
			new_line = new_line + line
	return new_line


	
def sorting(coords):
    new_value = ''
    myd = myDict()
    for line in coords.splitlines():
        element = line.split()
        element[0] = int(element[0])
        myd.add(element[0],line)
    for key, value in sorted(myd.iteritems()):
        value = value + '\n'
        new_value = new_value + value
    return new_value



def remove_line_numbers(sorted_coords):
    new_line = ''
    for line in sorted_coords.splitlines():
        element = line.split()
        line = line.replace(element[0], '')
        line = line + '\n'
        new_line = new_line + line
    return new_line



def add_atomic_numbers(pre_final_coords):
    count = 0
    new_line = ''
    for line in pre_final_coords.splitlines():
        if 'activation:' in line or 'lanczos:' in line:
            line = line.replace('activation:', '8' + '\n' + 'activation')
            line = line.replace('lanczos:', '8' + '\n' + 'lanczos')
        
        if not 'activation' in line and not 'lanczos' in line:
            count = count + 1
            if count==1 or count==5:
                line = '6' + '\t' + line
            else:
                line = '1' + '\t' + line
        if 'activation' in line or 'lanczos' in line:
            count = 0

        line = line + '\n'

        new_line = new_line + line

    return new_line

		

if __name__ == "__main__":
        initial_activation = coords_before_K0()
	activation = grep_activation_relaxation("SADDLE","K converged")
	relaxation = grep_activation_relaxation("K converged","SADDLE")

        init = remove_junk_lines(initial_activation, 'activation')
        act =  remove_junk_lines(activation, 'activation')
	lanc = remove_junk_lines(relaxation, 'lanczos')
       
        coords = init + act + lanc
        sorted_coords = sorting(coords)
        pre_final_coords = remove_line_numbers(sorted_coords)
        final_coords = add_atomic_numbers(pre_final_coords)

        print(final_coords)

