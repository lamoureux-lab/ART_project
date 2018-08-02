
import sys
from os import makedirs
from os.path import isfile, join, exists
from shutil import rmtree


def query_yes_no(question, default="yes"):
    """
    From Python recipe http://code.activestate.com/recipes/577058/
    Asks a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


def create_directory(directory):
    """
    Creates a new directory if it does not already exist
    :param directory:
    :return: boolean indicating that a new directory was created (True) or already existed and was not created (False)
    """
    if not exists(directory):
        makedirs(directory)
        return True
    return False


def check_file_exists(filename):
    """
    Checks that a file exists, printing a message that it is missing if it does not
    :param directory:
    :return: boolean indicating that the file does or does not exist
    """
    if not isfile(filename):
        print 'Missing: ' + filename
        return False
    return True


#
def check_for_missing_files_reset(file_directory, filename_list, message):
    '''
    Checks that an existing structure directory has the correct files to continue running ART
    The option will be given to skip directory and submit files for the remaining gaussian.inp
    or to clear the directory and restart

    :param file_directory:
    :param filename_list:
    :return:
    '''

    question = message
    reset = False
    file_missing = False

    for filename in filename_list:
        if not check_file_exists(join(file_directory, filename)):
            file_missing = True

    if file_missing:
        print 'A critical file is missing from: ' + file_directory + '\n'
        erase = query_yes_no(question)
        if erase:
            rmtree(file_directory)
            if create_directory(file_directory):
                reset = True

    return reset