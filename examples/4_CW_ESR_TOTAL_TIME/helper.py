import matplotlib.pyplot as plt
import matplotlib
import numpy as np

fontsize = 10
labelsize = 10
dpi = 72.27
print_dpi = 1200
fig_width_inches = 3.4
fig_height_inches = 3.4

def get_dc_data( filename='SpectraESR.dat' ):
   dat = np.loadtxt( filename )
   return dat

def plot_2d( ax_obj, x, y, xfact=1, yfact=1, swapxy=False, linewidth=1.5, marker='', color='k', label='', xlim=None, ylim=None, name=None ):
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

   ax_obj.plot( xplot, yplot, linewidth=linewidth, marker=marker, label=name )

def plot_esr_compare( datafiles=['rough/SpectraESR.dat','fine/SpectraESR.dat'], savename='SpectraESR.png', names=['rough','fine'], freq_diff=None):
   
   xlabel = r'Frequency (GHz)'
   ylabel = r'DC Current (pA)'
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
   freq = {}
   current = {}
   for i, file in enumerate(datafiles):
      esr_data[file] = get_dc_data( filename=file ) 
      freq[file] = esr_data[file][:,0]
      current[file] = esr_data[file][:,1]
      xmin = freq[file][0] 
      xmax = freq[file][-1]
      if freq_diff:
         xdiff = freq_diff
      else:
         xdiff = freq[file][1] - freq[file][0]
      plot_2d( ax, freq[file], current[file], xfact=xfact, yfact=yfact, xlim=[xmin,xmax,xdiff], marker='.', color=colors[i], name=names[i] )
   ax.legend(frameon=False)
   plt.tight_layout()
   plt.savefig( savename, bbox_inches='tight', dpi=print_dpi)


   
