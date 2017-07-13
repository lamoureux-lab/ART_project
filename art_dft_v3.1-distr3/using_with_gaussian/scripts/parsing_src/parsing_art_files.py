def set_env_config(natoms):
    """
    Sets critical parameters in gaussian_art.sh which contains the general configuration for the ART environment

    :param natoms:
    :return:
    """
    global script_directory
    global refconfig_filename
    global gaussian_art_filename

    with open(join(script_directory, gaussian_art_filename)) as input:

        updated_text = ''
        for line in input:
            original_line = line
            updated_line = ''

            #Remove comments from the strings
            line = line.split("#")
            comment_set = False
            if len(line) > 1:
                comment_set = True

            #Checks that GAU is set for energy calculation and sets it if it is not
            if line[0].find('setenv ENERGY_CALC') != -1:
                updated_line = 'setenv ENERGY_CALC GAU      '

            #Sets the number of atoms
            elif line[0].find('setenv NATOMS') != -1:
                updated_line = 'setenv NATOMS     '+ str(natoms) + '         '

            #Sets reference configuration file that is updated throughout
            elif line[0].find('setenv REFCONFIG') != -1:
                updated_line = 'setenv REFCONFIG        ' + refconfig_filename + '             '

            # Puts back configuration comments
            if comment_set and updated_line:
                updated_line = updated_line + '#' + line[1]
            elif updated_line:
                updated_line = updated_line + '\n'

            # updates the configuration file string
            if updated_line:
                updated_text = updated_text + updated_line
            else:
                updated_text = updated_text + original_line

    # overwrites the configuration file with appropriate values from the gaussian input file
    config = open(join(script_directory, gaussian_art_filename), 'w+')
    config.write(updated_text)
    config.close()