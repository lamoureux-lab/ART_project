'''Script to extract the activation and lanczos steps from the ART output.log file.
   Running this will provide the user with a clear picture of how many activation 
   steps it took before the algorithm changed to Lanczos and how many Lanczos steps 
   it took before the configuratuon converged to a saddle point and whether or not it
   converged.
'''

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


def add_saddle_line():
    with open('output.log') as f:
        new_line = ''
        count = 0
        for line in f:
            count = count + 1
            line = str(count) + '\t' + line
            if 'SADDLE' in line:
                new_line = new_line + line
        return new_line


def locate_sad_failed_lines(sorted_coords, nat):
    failed_list = []
    for line_number, line in enumerate(sorted_coords.splitlines()):
        if 'FAILED' in line:
            failed_list.append(int(line_number)-(int(nat)+1))
    return failed_list

def locate_sad_converged_lines(sorted_coords, nat):
    converged_list = []
    for line_number, line in enumerate(sorted_coords.splitlines()):
        if 'CONVERGED' in line:
            converged_list.append(int(line_number)-(int(nat)+1))
    return converged_list


def demarcate_lanczos(sorted_coords,failed_sad_coords, converged_sad_coords):
    new_line = ''
    for line_number, line in enumerate(sorted_coords.splitlines()):
        for i in range(0,len(failed_sad_coords)):
            if line_number == failed_sad_coords[i]:
                line = line.replace('lanczos', 'saddle_failed')
        for i in range(0, len(converged_sad_coords)):
            if line_number == converged_sad_coords[i]:
                line = line.replace('lanczos', 'saddle_converged')
        new_line = new_line + line + '\n'
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
			new_line = new_line + line + '\n'
        return new_line

	
def sorting(coords):
    new_value = ''
    mydict = {}
    for line in coords.splitlines():
        element = line.split()
        element[0] = int(element[0])
        mydict[element[0]] = line
    for key in sorted(mydict.keys()):
        mydict[key] = mydict[key] + '\n'
        new_value = new_value + mydict[key]    
    return new_value



def remove_line_numbers(demarcated_sad_coords):
    new_line = ''
    for line in demarcated_sad_coords.splitlines():
        element = line.split()
        line = line.replace(element[0], '')
        if not 'SADDLE' in line:
            new_line = new_line + line + '\n'
    return new_line

def fetch_atomic_numbers():
    with open('output.log') as f:
        new_atomic_list = []
        for line in f:
            if 'typat' in line:
                atomic_list = line.split()
                break
        for i in range(1,len(atomic_list)):
            new_atomic_list.append(atomic_list[i])
        
        return new_atomic_list

def fetch_number_of_atoms():
    with open('output.log') as f:
        for line in f:
            if 'nat' in line:
                element = line.split()
                break
        number_of_atoms = element[1]
        return number_of_atoms


def add_atomic_numbers(pre_final_coords, atomic_list, nat):
    count = 0
    new_line = ''
    for line in pre_final_coords.splitlines():
        if 'activation' in line or 'lanczos' in line or 'saddle_converged' in line or 'saddle_failed' in line:
            line = line.replace('activation', nat + '\n' + 'activation')
            line = line.replace('lanczos', nat + '\n' + 'lanczos')
            line = line.replace('saddle_converged', nat + '\n' + 'saddle_converged')
            line = line.replace('saddle_failed', nat + '\n' + 'saddle_failed') 

        if not 'activation' in line and not 'lanczos' in line and not 'saddle_converged' in line and not 'saddle_failed' in line:
            line = atomic_list[count] + '\t' + line
            count = count + 1

        if 'activation' in line or 'lanczos' in line  or 'saddle_converged' in line or 'saddle_failed' in line:
            count = 0

        new_line = new_line + line + '\n'
    return new_line

		

if __name__ == "__main__":
        initial_activation = coords_before_K0()
	activation = grep_activation_relaxation("SADDLE","K converged")
	relaxation = grep_activation_relaxation("K converged","SADDLE")

        init = remove_junk_lines(initial_activation, 'activation')
        act =  remove_junk_lines(activation, 'activation')
	lanc = remove_junk_lines(relaxation, 'lanczos')
        sad = add_saddle_line()

        nat = fetch_number_of_atoms()       
        atomic_list = fetch_atomic_numbers()

        coords = init + act + lanc + sad
        sorted_coords = sorting(coords) 
        failed_sad_coords = locate_sad_failed_lines(sorted_coords, nat)
        converged_sad_coords = locate_sad_converged_lines(sorted_coords, nat)
        demarcated_sad_coords = demarcate_lanczos(sorted_coords,failed_sad_coords,converged_sad_coords)
        pre_final_coords = remove_line_numbers(demarcated_sad_coords)
       
        final_coords = add_atomic_numbers(pre_final_coords,atomic_list,nat)
        print final_coords
       

