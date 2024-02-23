import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from tepy.globalvars import *
from tepy.tools import *

def get_dm( rho, dt, total_dim ):
   
   num_dim=4
   
   real_indices = np.array([ [get_py_indx_rhofile( i+1, j+1, 'real', total_dim ) for j in range(num_dim)] for i in range(num_dim) ])
   imag_indices = np.array([ [get_py_indx_rhofile( i+1, j+1, 'imag', total_dim ) for j in range(num_dim)] for i in range(num_dim) ])
   
   rho_data = np.zeros([num_dim,num_dim],dtype=complex)
   for i in range(num_dim):
      for j in range(num_dim):
         rho_data[i,j] = rho[dt,real_indices[i,j]] + 1j*rho[dt,imag_indices[i,j]]
   
   return rho_data

def plot_dmelements( num_dim=4, total_dim=8, datafile='RHO.dat', time_diff=None, vlines=None ):

   xlabel = r'Time (ns)'
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

   for i in range(num_dim):
      
      ylabel = r'Re[$\rho_{ij}$]'
      fig, ax = plt.subplots()
      fig.set_size_inches( fig_width_inches, fig_height_inches)
      fig.set_dpi(dpi)
      ax.tick_params(direction='in')
      ax.set_xlabel(xlabel, fontsize=fontsize)
      ax.set_ylabel(ylabel, fontsize=fontsize)

      real_indices = [get_py_indx_rhofile( i+1, j+1, 'real', total_dim ) for j in range(num_dim)]
      imag_indices = [get_py_indx_rhofile( i+1, j+1, 'imag', total_dim ) for j in range(num_dim)]

      for j in range(num_dim):
         plot_2d( ax, time, rho[:,real_indices[j]], xfact=xfact, yfact=yfact, xlim=[xmin,xmax,xdiff], ylim=[-1,1,0.5], marker=None, linewidth=1.5, color=colors[j], name=r'$\rho_{%d%d}$' % (i+1,j+1))
      if vlines:
         for line in vlines:
            ax.axvline( x=line, color='gray', linestyle='--', linewidth=0.75 )       
      
      ax.legend(frameon=False)
      plt.tight_layout()
      plt.savefig( 'rho_row%d_real.png' % (i), bbox_inches='tight', dpi=print_dpi)
      plt.close()
      
      ylabel = r'Im[$\rho_{ij}$]' 
      fig, ax = plt.subplots()
      fig.set_size_inches( fig_width_inches, fig_height_inches)
      fig.set_dpi(dpi)
      ax.tick_params(direction='in')
      ax.set_xlabel(xlabel, fontsize=fontsize)
      ax.set_ylabel(ylabel, fontsize=fontsize)
      
      for j in range(num_dim):
         plot_2d( ax, time, rho[:,imag_indices[j]], xfact=xfact, yfact=yfact, xlim=[xmin,xmax,xdiff], ylim=[-1,1,0.5], marker=None, linewidth=1.5, color=colors[j], name=r'$\rho_{%d%d}$' % (i+1,j+1))
      if vlines:
         for line in vlines:
            ax.axvline( x=line, color='gray', linestyle='--', linewidth=0.75 ) 
           
      ax.legend(frameon=False)
      plt.tight_layout()
      plt.savefig( 'rho_row%d_imag.png' % (i), bbox_inches='tight', dpi=print_dpi)
      plt.close()

def plot_pops( num_populations=3, datafiles=['POPULATIONS.dat'], savename='Populations.png', names=[r'$\downarrow$',r'$\uparrow$',r'$\emptyset$'], time_diff=None, vlines=None ):

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
      if vlines:
         for line in vlines:
            ax.axvline( x=line, color='gray', linestyle='--', linewidth=0.75 )
   ax.legend(frameon=False)
   plt.tight_layout()
   plt.savefig( savename, bbox_inches='tight', dpi=print_dpi)
      
def check_final_states( filename, times):
   rho = get_data( filename=filename )
   arraysize = len(rho[0,:])
   print(len(rho[:,0]))
   diag_states = [(1,1),(2,2),(3,3),(4,4)]
   coh_states = [(1,4),(4,1)]
   other_states = []
   for istate in range(1,9):
      for jstate in range(1,9):
         state = (istate, jstate)
         if state not in diag_states:
            if state not in coh_states:
               other_states.append((istate,jstate))
            
   for time in times:
      print("Time (ns): ",time)
      itime = int(time / (1000./500000))-1
      diag_val = 0.0
      coh_val = 0.0
      other_val = 0.0
      for state in diag_states:
         istate, jstate = state
         i = get_py_indx_rhofile( istate, jstate, 'real', 8 )
         j = get_py_indx_rhofile( istate, jstate, 'imag', 8 )
         diag_val = diag_val + abs(rho[itime,i] + (1j*rho[itime,j]**2.))
      for state in coh_states:
         istate, jstate = state
         i = get_py_indx_rhofile( istate, jstate, 'real', 8 )
         j = get_py_indx_rhofile( istate, jstate, 'imag', 8 )
         coh_val = coh_val + abs(rho[itime,i] + (1j*rho[itime,j]**2.))
      for state in other_states:
         istate, jstate = state
         i = get_py_indx_rhofile( istate, jstate, 'real', 8 )
         j = get_py_indx_rhofile( istate, jstate, 'imag', 8 )
         other_val = other_val + abs(rho[itime,i] + (1j*rho[itime,j]**2.))
      print('Diag vals:',diag_val)
      print('Coh vals:',coh_val)
      print('Other vals:',other_val)
      #if abs(other_val) >= 0.1:
      #print('Other States %d %d:' % (istate,jstate),val)