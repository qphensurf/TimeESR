from helper import *
from os import sys

plot_esr_compare( freq_diff=0.1, datafiles=[sys.argv[1],sys.argv[2]], names=[sys.argv[3],sys.argv[4]], savename=sys.argv[5] )
