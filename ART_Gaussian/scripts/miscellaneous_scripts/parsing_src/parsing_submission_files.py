
from os.path import join, dirname

def set_submission_script(script_directory, submission_script, structure_output_directory, link0_section,
                          optimize=False, time = None, memory = None, processor_info = None):
    """
        Sets submissions parameters in gaussian_art.sh which contains the general configuration for the ART environment

        :param natoms:
        :return:
        """
    # global script_directory
    # submission_script_template = join(dirname(__file__), submission_script)
    submission_script_template = join(dirname(__file__), join(script_directory, submission_script))
    link0_arr = link0_section.splitlines()

    with open(submission_script_template) as input:
        updated_text = ''
        for line in input:
            original_line = line
            updated_line = ''

            # Sets the time
            time_flag = '#PBS -l walltime'
            if time_flag in line:
                if time:
                    updated_line = '#PBS -l ' + time + '\n'

            # Sets the memory
            memory_flag = '#PBS -l mem'
            if memory_flag in line:
                if memory:
                    updated_line = '#PBS -l ' + memory + '\n'
                elif optimize:
                    for item in link0_arr:
                        link0_memory = extract_link0_parameter('mem', item)
                        #Checks if the memory is set in the gaussian input file
                        if link0_memory:
                            updated_line = '#PBS -l mem=' + link0_memory + '\n'

            # Sets the processor information
            processor_flag = '#PBS -l nodes'
            if processor_flag in line:
                if processor_info:
                    updated_line = '#PBS -l ' + processor_info + '\n'
                elif optimize:
                    # % nproc = 12
                    for item in link0_arr:
                        link0_nproc = extract_link0_parameter('nproc', item)
                        # Checks if the nproc is set in the gaussian input file
                        if link0_nproc:
                            updated_line = '#PBS -l nodes=1:ppn=' + link0_nproc + '\n'

            # updates the submission file string
            if updated_line:
                updated_text = updated_text + updated_line
            else:
                updated_text = updated_text + original_line

    # overwrites the submissions file with appropriate values from the gaussian input file
    sub = open(join(structure_output_directory, submission_script), 'w+')
    sub.write(updated_text)
    # print updated_text
    sub.close()

def extract_link0_parameter(flag, link0_line):
    """
    Gets the parameter assigned for memory or number of processors in a gaussian.inp file
    :return:
    """
    if flag in link0_line:
        return link0_line.split("=")[1].strip()
    return None
