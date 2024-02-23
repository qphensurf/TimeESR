from os import sys
from tepy.density import *

numstates = int(sys.argv[1])
filename = sys.argv[2]
total_dim = int(sys.argv[3])
time_diff = int(sys.argv[4])
if sys.argv[5]:
   line_nums = int(sys.argv[5])
   lines = sys.argv[6:6+line_nums]
   vlines = [float(line) for line in lines]
plot_dmelements( time_diff=time_diff, num_dim=numstates, datafile=filename, total_dim=total_dim, vlines=vlines )
