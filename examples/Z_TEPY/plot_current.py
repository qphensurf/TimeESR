from os import sys
from tepy.current import *

filename = sys.argv[1]
outputname = sys.argv[2]
time_diff = int(sys.argv[3])
if sys.argv[4]:
   line_nums = int(sys.argv[4])
   lines = sys.argv[5:5+line_nums]
   vlines = [float(line) for line in lines]
plot_current( time_diff=time_diff, datafile=filename, savename=outputname, vlines=vlines )
