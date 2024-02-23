import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from tepy.globalvars import *
from tepy.tools import *

def plot_current( datafile='C.dat', savename='Current.png', time_diff=None, vlines=None ):

   xlabel = r'Time (ns)'
   ylabel = r'Current (pA)'
   xfact = 1
   yfact = 1
   
   matplotlib.rc('xtick', labelsize=labelsize)
   matplotlib.rc('ytick', labelsize=labelsize)
   plt.rcParams['figure.figsize'] = (fig_width_inches, fig_height_inches)
   plt.rcParams['figure.dpi'] = dpi

   prop_cycle = plt.rcParams['axes.prop_cycle']
   colors = prop_cycle.by_key()['color']

   fig, ax = plt.subplots()
   fig.set_size_inches( fig_width_inches, fig_height_inches)
   fig.set_dpi(dpi)
   ax.tick_params(direction='in')
   ax.set_xlabel(xlabel, fontsize=fontsize)
   ax.set_ylabel(ylabel, fontsize=fontsize)
   
   esr_data = {}
   time = {}
   file = datafile
   esr_data[file] = get_data( filename=file )
   time[file] = esr_data[file][:,0]
   xmin = time[file][0] 
   xmax = time[file][-1]
   if time_diff:
      xdiff = time_diff
   else:
      xdiff = time[file][1] - time[file][0]
   plot_2d( ax, time[file], esr_data[file][:,1], xfact=xfact, yfact=yfact, xlim=[xmin,xmax,xdiff], marker=None, linewidth=1.5, color=colors[0] )
   if vlines:
      for line in vlines:
         ax.axvline( x=line, color='gray', linestyle='--', linewidth=0.75 )
   plt.tight_layout()
   plt.savefig( savename, bbox_inches='tight', dpi=print_dpi)