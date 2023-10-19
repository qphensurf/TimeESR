from helper import *
from os import sys

numstates = int(sys.argv[1])
filename = sys.argv[2]
outputname = sys.argv[3]
names = sys.argv[4:4+numstates]
plot_spinpops( time_diff=30, num_populations=numstates, datafiles=[filename], names=names, savename=outputname )
