import re
from os.path import join

class load_input:
    """
    Loads the values gaussian parameters into a dictionary object
    :param gaussian_input_params
    :return:
    """
    def __init__(self, input_file):
        self.input = input_file
        self.gaussian_input_params = {'link0_section': '', 'route_section': '',
                                      'title': '', 'natoms': 0, 'charge': None, 'multiplicity': None, 'atom_coordinates' : ''}
        self.gaussian_input_params = self.read_gaussian_input()

    def read_gaussian_input(self):
        input_file = self.input
        params = self.gaussian_input_params
        with open(input_file) as input:

            #The periodic_table_dict will be loaded once if conversion is needed from atomic symbol to number

            section_number = 1
            for line in input:
                    # This will load the relevant information from a gaussian input file expecting:
                    # 1. A number of % lines (or maybe none) representing the link 0 section
                    # 2. A number of # line describing the job representing the route section
                    # 3. Blank line, title line, blank line
                    # 4. Charge + Spin multiplicity line (two integers)
                    # 5. Atom coordinates

                    # Remove all comments from the gaussian input file
                    comment_set = False
                    line = line.split("!")
                    if len(line) > 1:
                        comment_set = True
                    line = line[0]
                    if comment_set:
                        line = line + '\n'

                    #Check for an empty line to increment section
                    if line.strip() == '':
                        if not comment_set:
                            section_number += 1
                        continue

                    #Link0 and route section
                    if section_number == 1:
                        #For case issues
                        line = line.lower()
                        #Link0
                        if line.find('%') != -1:
                            params['link0_section'] = params['link0_section'] + line + " "
                        #Route section
                        elif line.find('#') != -1:
                            params['route_section'] = params['route_section'] + line

                    #Title section
                    elif section_number == 2:
                        params['title'] = line.strip()

                    #Molecule specification section
                    if section_number == 3:
                        string_integers = line.split()
                        params['charge'] = int(string_integers[0])
                        params['multiplicity'] = int(string_integers[1])
                        section_number += 1
                        continue

                    #Atom coordinates section
                    if section_number == 4:
                        params['natoms'] = params['natoms'] + 1
                        params['atom_coordinates'] = params['atom_coordinates'] + line

                    if section_number == 5:
                        break
            return params

    def symbol_to_atomic_number(self, periodic_table):
        # Checks if the atom is in non-numerical format (e.g., H instead of 1)
        periodic_table_dict = None
        atom_coordinates = self.gaussian_input_params['atom_coordinates'].strip().split('\n')
        self.gaussian_input_params['atom_coordinates'] = ''


        for line in atom_coordinates:
            temp_line = line.split()
            element = temp_line[0]
            if not element.isdigit():
                # Loads the periodic table data if it is not already set
                if not periodic_table_dict:
                    periodic_table_dict = self._load_periodic_table(periodic_table)
                # Replaces atomic symbols with atomic numbers
                temp_line[0] = periodic_table_dict[element]
                line = '    '.join(temp_line) + '\n'
            self.gaussian_input_params['atom_coordinates'] = self.gaussian_input_params['atom_coordinates'] + line

    def _load_periodic_table(self, periodic_table):
        """
        Loads a periodic table csv file to extract the atomic symbols and numbers to create a mapping between these
        :return: A dictionary mapping atomic symbol to number
        """

        atomic_symbol_number_map = {}
        with open(periodic_table) as input:
            for line in input:
                line = line.split(',')
                atomic_symbol_number_map[line[1].strip()] = line[0].strip()

        return atomic_symbol_number_map

    def remove_route_parameter(self, parameter_list):
        # Removes options such as opt and force, which will be put back in by the ART code at different stages

        route_section = self.gaussian_input_params['route_section']
        self.gaussian_input_params['route_section'] = ''

        for line in route_section:
            for option in parameter_list:
                line = self._option_removal_helper(option, line)
            if line.strip() != '#':
                self.gaussian_input_params['route_section'] = self.gaussian_input_params['route_section'] + line


    def _option_removal_helper(self, option_type, line):
            # Handles different options to keyword structures for removal to have them dynamically assigned by ART:
            # keyword = option
            # keyword(option)
            # keyword=(option1, option2)
            # keyword(option1, option2)

            # For keyword = option
            line = re.sub(r'(\s)' + re.escape(option_type) + '=\w+', r'\1', line)
            line = re.sub(r'(#)' + re.escape(option_type) + '=\w+', r'\1', line)
            line = re.sub(r'(\s)' + re.escape(option_type) + ' =\w+', r'\1', line)
            line = re.sub(r'(#)' + re.escape(option_type) + ' =\w+', r'\1', line)
            line = re.sub(r'(\s)' + re.escape(option_type) + ' = \w+', r'\1', line)
            line = re.sub(r'(#)' + re.escape(option_type) + ' = \w+', r'\1', line)
            line = re.sub(r'(\s)' + re.escape(option_type) + ' = \w+', r'\1', line)
            line = re.sub(r'(#)' + re.escape(option_type) + ' = \w+', r'\1', line)

            # For keyword(option..)
            line = re.sub(r'(\s)' + re.escape(option_type) + '\([^()]*\)', '', line)
            line = re.sub(r'(#)' + re.escape(option_type) + '\([^()]*\)', '', line)

            # For keyword=(option..)
            line = re.sub(r'(\s)' + re.escape(option_type) + '=\([^()]*\)', '', line)
            line = re.sub(r'(#)' + re.escape(option_type) + '=\([^()]*\)', '', line)

            # For lone keyword
            line = re.sub(r'(\s)' + re.escape(option_type), r'\1', line)
            line = re.sub(r'(#)' + re.escape(option_type), r'\1', line)

            return line

    def extract_link0_parameter(flag, link0_line):
        """
        Gets the parameter assigned for memory or number of processors in a gaussian.inp file
        :return:
        """
        if flag in link0_line:
            return link0_line.split("=")[1].strip()
        return None