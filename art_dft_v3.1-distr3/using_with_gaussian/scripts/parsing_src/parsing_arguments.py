
import argparse

def argument_parser():
    parser = argparse.ArgumentParser(description='Gaussian_ART ... ') #TODO add a good application description
    parser.add_argument('-d', '--project_directory',
                        help='project directory containing gaussian input files', default=default_input_file_directory)
    parser.add_argument('-f', '--input_files', nargs='*',
                        help='specific input files to submit from project directory')
    parser.add_argument('-s', '--submission_type', choices=['GREX', 'PSI'], default=default_submission,
                        help='GREX or PSI submissions supported. ' + default_submission + ' is the default')
    parser.add_argument('-r', '--reset', action='store_true',
                        help='resets the output directory instead of continuing from latest referenced configuration')
    parser.add_argument('-t', '--time',
                        help='changes the submissions job time (e.g., -t \'walltime=010:00:00\' ')
    parser.add_argument('-m', '--memory',
                        help='changes the submissions job memory cap (e.g., -m \'mem=8000MB\' ')
    parser.add_argument('-p', '--processor_info',
                        help='changes the submissions job processor cap (e.g., -p \'nodes=1:ppn=8\' ')
    parser.add_argument('-o', '--optimize', action='store_true',
                        help='sets submission job memory and number of processors to what is specified in the file\n '
                             '(note: specifying --memory or --processor_info will take precedence)')
    # args parser.parse_args()
    return parser
