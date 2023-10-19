from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from helper import *

# ----- START User options ----- #

# Number of steps to skip in the data file for each picture added to the movie
slice_step = 500

# Number of printed states in RHO.dat
rho_total_printed_states = 3

# Chosen eigenstates for the Bloch vector (start counting from 1)
rho_0_state_num = 1
rho_1_state_num = 2

# Rotating frame transformation
rotframe = True  # Enable transformation
f_d = 13.995322  # Driving frequency, including de-tuning (GHz)
eps1 = 0.0       # Energy of level 1 (GHz)
eps2 = 13.995322 # Energy of level 2 (GHz)

# Frequency analysis
print_freqs = True

# ----- END User options ----- #

# Load projected Bloch sphere data file (0: time, 1:4: vector x, y, z)
data = project_bloch_sphere( rho_0_state_num, rho_1_state_num, rho_total_printed_states, rotframe, f_d, eps1, eps2, './RHO.dat', print_freqs )

# Save picture of Bloch vector trajectory and animate Bloch vector over entire trajectory
animate_bloch_plot( data[:,1:4], data[:,0], slice_step )

