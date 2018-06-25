import argparse

parser = argparse.ArgumentParser(description = 'ART coupled with Gaussian')

parser.add_argument('-k', '--keyword', help = 'Route section keyword')

args = parser.parse_args()


keyword = args.keyword

input_file =                 # Filled in when the artgau.py is called

with open(input_file) as f:
    with open('art2gaussian.inp','w+') as m:
        i = 0
        for index, line in enumerate(f,1):
            if '%' in line or '#' in line:
                i = i + 1
            if '#' in line:
                if keyword == 'opt':
                    line = line.replace('force', 'opt')
            if (index <= (i + 4)):
                m.write(line)
