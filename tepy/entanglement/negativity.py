import numpy as np
from scipy.linalg import eigvals
from tepy.globalvars import *
from tepy.tools import *
from tepy.density import *
  
def negativity( rho ):
   """
   Negativity of subsystem rho
   N(rho) = \abs( \sum_lambda_i<0 labmda_i )
   = \frac{1}{2} \sum_i
   \abs{labmda_i} - \lambda_i,  
   where lambda_i are the eigenvalues of 
   \rho^\Gamma_A, which is the partial transpose of \rho with respect to subsystem A
   See https://toqito.readthedocs.io/en/latest/_autosummary/toqito.channels.partial_transpose.html for an example
   """
   row_0 = [rho[0,0],rho[0,1],rho[2,0],rho[2,1]]
   row_1 = [rho[1,0],rho[1,1],rho[3,0],rho[3,1]]
   row_2 = [rho[0,2],rho[0,3],rho[2,2],rho[2,3]]
   row_3 = [rho[1,2],rho[1,3],rho[3,2],rho[3,3]]
   rho_gamma_A = np.array([row_0,row_1,row_2,row_3],dtype=complex)
   evals = eigvals(rho_gamma_A).real
   N = ((abs(evals)-evals)/2).sum()
   return N

def plot_negativity( total_dim=8, datafile='RHO.dat', savename='negativity.png', time_diff=None, vlines=None ):
   
   xlabel = r'Time (ns)'
   ylabel = r'Negativity'
   xfact = 1
   yfact = 1
      
   matplotlib.rc('xtick', labelsize=labelsize)
   matplotlib.rc('ytick', labelsize=labelsize)
   plt.rcParams['figure.figsize'] = (fig_width_inches, fig_height_inches)
   plt.rcParams['figure.dpi'] = dpi

   prop_cycle = plt.rcParams['axes.prop_cycle']
   colors = prop_cycle.by_key()['color']
   
   rho = get_data( filename=datafile )
   time = rho[:,0]
   xmin = time[0] 
   xmax = time[-1]
   if time_diff:
      xdiff = time_diff
   else:
      xdiff = time[1] - time[0]

   fig, ax = plt.subplots()
   fig.set_size_inches( fig_width_inches, fig_height_inches)
   fig.set_dpi(dpi)
   ax.tick_params(direction='in')
   ax.set_xlabel(xlabel, fontsize=fontsize)
   ax.set_ylabel(ylabel, fontsize=fontsize)
   
   negat = np.zeros(len(time))
   for dt in range(len(time)):
      rho_data = get_dm( rho, dt, total_dim )
      negat[dt] = negativity( rho_data )
            
   plot_2d( ax, time, negat, xfact=xfact, yfact=yfact, xlim=[xmin,xmax,xdiff], ylim=[0,0.5,0.1], marker=None, linewidth=1.5, color=colors[0])
   
   if vlines:
      for line in vlines:
         ax.axvline( x=line, color='gray', linestyle='--', linewidth=0.75 )       
      
   plt.tight_layout()
   plt.savefig( savename, bbox_inches='tight', dpi=print_dpi)  