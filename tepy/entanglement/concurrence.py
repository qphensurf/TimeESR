import numpy as np
from scipy.linalg import sqrtm, eigvals
from tepy.globalvars import *
from tepy.tools import *
from tepy.density import *

def concurrence( rho ):
   """
   Concurrence of density matrix rho
   C(rho) = max(0,lambda_1-lambda_2-lambda_3-lambda_4)
   where lambda_i are the eigenvalues of 
   R = sqrt( sqrt(rho) rho_tilde sqrt(rho) )
   where rho_tilde = (sy tensor sy) rho^* (sy tensor sy)
   """
   sysy = np.array([[0.,0.,0.,-1.],[0.,0.,1.,0.],[0.,1.,0.,0.],[-1.,0.,0.,0.]],dtype=complex)
   rho_tilde = np.matmul( sysy, np.matmul( rho.conj(), sysy ) )
   rho_sqrt = sqrtm(rho)
   R = sqrtm( np.matmul( rho_sqrt, np.matmul( rho_tilde, rho_sqrt ) ) )
   evals = eigvals(R).real
   evals_sorted = np.sort(evals)[::-1]
   return np.max([0.,evals_sorted[0]-evals_sorted[1]-evals_sorted[2]-evals_sorted[3]])
   
def plot_concurrence( total_dim=8, datafile='RHO.dat', savename='concurrence.png', time_diff=None, vlines=None ):
   
   xlabel = r'Time (ns)'
   ylabel = r'Concurrence'
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
   
   concur = np.zeros(len(time))
   for dt in range(len(time)):
      rho_data = get_dm( rho, dt, total_dim )
      concur[dt] = concurrence( rho_data )
            
   plot_2d( ax, time, concur, xfact=xfact, yfact=yfact, xlim=[xmin,xmax,xdiff], ylim=[0,1,0.2], marker=None, linewidth=1.5, color=colors[0])
   
   if vlines:
      for line in vlines:
         ax.axvline( x=line, color='gray', linestyle='--', linewidth=0.75 )       
      
   plt.tight_layout()
   plt.savefig( savename, bbox_inches='tight', dpi=print_dpi)  