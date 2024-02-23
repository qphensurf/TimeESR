from os import sys
from tepy.negativity import *

filename = sys.argv[1]
total_dim = int(sys.argv[2])
outputname = sys.argv[3]
time_diff = int(sys.argv[4])
if sys.argv[5]:
   line_nums = int(sys.argv[5])
   lines = sys.argv[6:6+line_nums]
   vlines = [float(line) for line in lines]
plot_negativity( total_dim=total_dim, datafile=filename, savename=outputname, time_diff=time_diff, vlines=vlines )