

input_file = 'art2gaussian.inp'

def create_ref_config():
    #TODO find out how to get:
    # total_energy:   -2144.09236620000
    #S   100.000000000000        100.000000000000        100.000000000000
    #Is this obtained from initial gaussian run?
    pass

def set_env_config():
    #TODO
    pass

def load_input(gaussian_input_params):

    params = gaussian_input_params
    counter = 0
    section_number = 1

    with open(input_file) as input:

        for line in input:
                # This will load the relevant information from a gaussian input file expecting:
                # 1. A number of % lines (or maybe none)
                # 2. One # line describing the job
                # 3. Blank line, title line, blank line
                # 4. Charge + Spin multiplicity line (two integers)
                # 5. Atom coordinates

                #Check for an empty line to increment section
                print line
                if line.strip() == '':
                    section_number += 1
                    continue

                #Check whether the line is in the first section
                #TODO find out whether to expect more than these two % lines for other_details section
                if section_number == 1:
                    if line.find('%mem') != -1:
                        params['mem'] = line
                    elif line.find('%nproc') != -1:
                        params['nproc'] = line
                    elif line.find('%') != -1:
                        params['other_details'] = params['other_details'] + line + " "
                    elif line.find('#') != -1:
                        params['description'] = line

                #Title section
                elif section_number == 2:
                    params['title'] = line

                #Charge + Spin multiplicity section
                if section_number == 3:
                    string_integers = line.split()
                    params['charge'] = int(string_integers[0])
                    params['multiplicity'] = int(string_integers[1])
                    section_number += 1

                #Atom coordinates
                if section_number == 4:
                    params['natoms'] = params['natoms'] + 1
                    params['atom_coordinates'] = params['atom_coordinates'] + line

                if section_number == 5:
                    break

        return params



if __name__ == "__main__":


    gaussian_input_params = {'mem': '', 'nproc': '', 'other_details': '', 'description': '',
                    'title': '', 'natoms': 0, 'charge': None, 'multiplicity': None, 'atom_coordinates' : ''}
    print load_input(gaussian_input_params)
    # print 'nproc: ' + nproc
    # print