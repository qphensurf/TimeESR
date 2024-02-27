from os import sys
from tepy.density import *

numstates = int(sys.argv[1])
filename = sys.argv[2]
total_dim = int(sys.argv[3])
time_diff = int(sys.argv[4])
rot_freq = float(sys.argv[5])
freq_width = float(sys.argv[6])
nfreq = int(sys.argv[7])
state1 = int(sys.argv[8])
state2 = int(sys.argv[9])
if sys.argv[10]:
   line_nums = int(sys.argv[10])
   lines = sys.argv[11:11+line_nums]
   vlines = [float(line) for line in lines]
plot_dmelements( time_diff=time_diff, num_dim=numstates, datafile=filename, 
                 total_dim=total_dim, vlines=vlines, rot=True, rot_freq=rot_freq, 
                 freq_width=freq_width, nfreq=nfreq, state1=state1, state2=state2 )