from os import sys
from tepy.density import *

numstates = int(sys.argv[1])
filename = sys.argv[2]
outputname = sys.argv[3]
names = sys.argv[4:4+numstates]
time_diff = int(sys.argv[4+numstates])
if sys.argv[4+numstates+1]:
   line_nums = int(sys.argv[4+numstates+1])
   lines = sys.argv[4+numstates+2:4+numstates+2+line_nums]
   vlines = [float(line) for line in lines]
plot_pops( time_diff=time_diff, num_populations=numstates, datafiles=[filename], names=names, savename=outputname, vlines=vlines )
