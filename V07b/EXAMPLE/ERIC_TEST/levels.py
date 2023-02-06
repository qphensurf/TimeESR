import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy import integrate
import sys
import os

# All values in meV
left_elec = -1.0
right_elec = 1.0
levels = [ 0.000, 0.05788, 0.07856, 0.1158, 10.03, 10.09, 90.03, 90.09 ]

# Get levels
#try:
#   os.popen("grep -a 'number of electrons' "+sys.argv[1] ).read()
#   nelec=float(os.popen("grep -a 'number of electrons' "+sys.argv[1] ).read().split()[4])
#   print("Number of electrons = ", nelec)
#except:
#   print("CRITICAL ERROR: Number of electrons not found!")
#   exit()

#nelecint = int(nelec)
#band_plot_limit_lower = -1.0
#band_plot_limit_upper = 0.5
plot_bands = True

glb_fontsize = 10
glb_labelsize = 10
glb_dpi = 72.27
print_dpi = 1200 
glb_figwdth_inches = 4.0
glb_fighght_inches = 4.0
glb_xaxis_unit = ''
glb_yaxis_unit = ' (meV)'

# set figure defaults
matplotlib.rc('xtick', labelsize=glb_labelsize)
matplotlib.rc('ytick', labelsize=glb_labelsize)
plt.rcParams["figure.figsize"] = (glb_figwdth_inches,glb_fighght_inches)
plt.rcParams["figure.dpi"] = glb_dpi

# setup colors
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

# plot
if plot_bands:
   fig1, ax1 = plt.subplots()
   fig1.set_size_inches(glb_figwdth_inches,glb_fighght_inches)
   fig1.set_dpi(glb_dpi)
   ax1.margins(x=0)
   ax1.set_ylabel(r'$E$'+glb_xaxis_unit, fontsize=glb_fontsize)
   ax1.tick_params(direction='in')
   for i in range( len( levels) ):
       ax1.plot( [1, 5], [levels[i], levels[i]], '-', lw = 1, color="C0" )
   ax1.set_xticks([],[])
   #ax1.set_ylim( band_plot_limit_lower, band_plot_limit_upper )
   fig1.tight_layout()
   plt.figure(fig1)
   plt.savefig('./LEVELS.png',bbox_inches="tight",dpi=print_dpi)
