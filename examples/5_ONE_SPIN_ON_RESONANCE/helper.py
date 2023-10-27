import matplotlib.pyplot as plt
import matplotlib
import numpy as np

fontsize = 10
labelsize = 10
dpi = 72.27
print_dpi = 1200
fig_width_inches = 3.4
fig_height_inches = 3.4

def get_data( filename='POPULATIONS.dat' ):
   dat = np.loadtxt( filename )
   return dat

def plot_2d( ax_obj, x, y, xfact=1, yfact=1, swapxy=False, linewidth=1.5, marker='', color='k', label='', xlim=None, ylim=None, name=None, markersize=6 ):
   # limits are in the form [min, max, diff]   

   if swapxy:
      xplot, yplot = y * float(yfact), x * float(xfact)
   else:
      xplot, yplot = x * float(xfact), y * float(yfact)
   
   if xlim:
      ticks = np.arange(xlim[0],xlim[1]+(xlim[2]/100),xlim[2])
      ax_obj.set_xlim(xlim[0],xlim[1])
      ax_obj.set_xticks(ticks)

   if ylim:
      ticks = np.arange(ylim[0],ylim[1]+(ylim[2]/100),ylim[2])
      ax_obj.set_ylim(ylim[0],ylim[1])
      ax_obj.set_yticks(ticks)

   ax_obj.plot( xplot, yplot, linewidth=linewidth, marker=marker, label=name, markersize=markersize )

def plot_spinpops( num_populations=3, datafiles=['POPULATIONS.dat'], savename='SpinPopulations.png', names=[r'$\downarrow$',r'$\uparrow$',r'$\emptyset$'], time_diff=None ):

   xlabel = r'Time (ns)'
   ylabel = r'Population'
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
   for i, file in enumerate(datafiles):
      esr_data[file] = get_data( filename=file )
      time[file] = esr_data[file][:,0]
      xmin = time[file][0] 
      xmax = time[file][-1]
      if time_diff:
         xdiff = time_diff
      else:
         xdiff = time[file][1] - time[file][0]
      for j in range(num_populations):
         plot_2d( ax, time[file], esr_data[file][:,1+j], xfact=xfact, yfact=yfact, xlim=[xmin,xmax,xdiff], marker=None, linewidth=1.5, color=colors[j], name=names[j] )
   ax.legend(frameon=False)
   plt.tight_layout()
   plt.savefig( savename, bbox_inches='tight', dpi=print_dpi)


   
