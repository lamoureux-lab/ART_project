from subprocess import call
import glob
files = glob.glob(min1*)
for file in files:
	call(['qsub', '-N','gau_art_'+ structure, grex_submission_script], shell=False)
