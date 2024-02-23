import numpy as np
from scipy.linalg import sqrtm, eigvals
from tepy.globalvars import *
from tepy.tools import *
from tepy.density import *

def fidelity( rho_1, rho_2 ):
   """
   Fidelity of two density matrices rho_1 and rho_2
   F = ( Tr [ Sqrt[ Sqrt[ dm_1 ] . dm_2 . Sqrt[ dm_1 ] ] ] )^2
   F = ( Tr [ Sqrt[ dm_1 . dm_2 ] ] ^2
   F should be between 0 and 1, with 1 indicating that the two density matrices are the same
   """
   return (np.trace(sqrtm(np.matmul(rho_1,rho_2))))**2.0
   
def fidelity_qutip( rho_a, rho_b, pure_a=False, pure_b=False ):
   """
   Fidelity of two density matrices rho_1 and rho_2 using the method in Qutip to avoid NaNs and better stabilize the sqrt of a matrix
   Main differences: swap the order of pure state B so that we can take a more numerically
            # stable square root of B.
   """   
   
   if pure_a:
      sqrtmA = rho_a
   else:
      if pure_b:
         rho_b, rho_a = rho_a, rho_b
      sqrtmA = sqrtm(rho_a)
   
   eig_vals = eigvals(sqrtmA * rho_b * sqrtmA)
   return float(np.real(np.sqrt(eig_vals[eig_vals > 0]).sum()))
   
   return 
   
def plot_fidelity( total_dim=8, datafile='RHO.dat', savename='fidelity.png', time_diff=None, vlines=None, belltype='phiplus', qutip_method=False ):

   xlabel = r'Time (ns)'
   ylabel = r'Fidelity'
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
   
   if belltype == 'phiplus':
      # 00 + 11
      rho_exact = 0.5*np.array([[1.,0.,0.,1.],[0.,0.,0.,0.],[0.,0.,0.,0.],[1.,0.,0.,1.]],dtype=complex)
   elif belltype == 'phiminus':
      # 00 - 11
      rho_exact = 0.5*np.array([[1.,0.,0.,-1.],[0.,0.,0.,0.],[0.,0.,0.,0.],[-1.,0.,0.,1.]],dtype=complex)
   elif belltype == 'psiplus':
      # 01 + 10
      rho_exact = 0.5*np.array([[0.,0.,0.,0.],[0.,1.,1.,0.],[0.,1.,1.,0.],[0.,0.,0.,0.]],dtype=complex)
   elif belltype == 'psiminus':
      # 01 - 10
      rho_exact = 0.5*np.array([[0.,0.,0.,0.],[0.,1.,-1.,0.],[0.,-1.,1.,0.],[0.,0.,0.,0.]],dtype=complex)
      
   fid = np.zeros(len(time))
   for dt in range(len(time)):
      rho_data = get_dm( rho, dt, total_dim )
      if qutip_method:
         fid[dt] = np.real(fidelity_qutip(rho_exact, rho_data, pure_a=False, pure_b=True))
      else:
         fid[dt] = np.real(fidelity(rho_exact, rho_data))
      #if np.isnan(fid[dt]):
      #   print("Possible non-Hermitian rho_data found for dt =",dt)
   
   plot_2d( ax, time, fid, xfact=xfact, yfact=yfact, xlim=[xmin,xmax,xdiff], ylim=[0,1,0.2], marker=None, linewidth=1.5, color=colors[0])

   '''
   x_curve, y_curve, tau = get_neg_exp_curve( time, fid )
   
   print("decay time: ",tau)
   
   ax.plot( x_curve, y_curve, linewidth=1.5, color='black', linestyle = '--', marker=None )
   '''
   if vlines:
      for line in vlines:
         ax.axvline( x=line, color='gray', linestyle='--', linewidth=0.75 )       
      
   plt.tight_layout()
   plt.savefig( savename, bbox_inches='tight', dpi=print_dpi) 