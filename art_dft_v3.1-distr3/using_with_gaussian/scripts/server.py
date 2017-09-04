import os
import sys
import subprocess
import time

if __name__=='__main__':
    for i in range(0,4):
        time.sleep(3)
        p = subprocess.Popen("python worker.py %d"%i,shell=True)