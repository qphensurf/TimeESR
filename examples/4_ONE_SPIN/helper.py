import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.fft import fft, fftfreq
from qutip import *
from mpl_toolkits.mplot3d import Axes3D

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()
        
def extract_peak_frequency(data, start, end, steps):
    sample_rate = (end-start)/steps
    fftfrequencies = fftfreq(steps,sample_rate)[:steps//2]
    fft_data = (2.0/steps) * np.abs(fft(data)[:steps//2])
    peak_coefficient = np.argmax(np.abs(fft_data))
    peak_freq = fftfrequencies[peak_coefficient]
    
    return abs(peak_freq)
    
def get_fft_spectra(data, start, end, steps):
    sample_rate = (end-start)/steps
    fftfrequencies = fftfreq(steps,sample_rate)[:steps//2]
    fft_data = (2.0/steps) * np.abs(fft(data)[:steps//2])
    
    return fft_data, fftfrequencies

def get_py_indx_rhofile( istate, jstate, imag_or_real, total_states ):
   # State number should be given in 1 base

   # Example for 3 state system
   # python index  0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18
   # column type   t Re 11 Re 21 Re 31 Re 12 Re 22 Re 32 Re 13 Re 23 Re 33 Im 11 Im 21 Im 31 Im 12 Im 22 Im 32 Im 13 Im 23 Im 33

   # Converted to python index format (note the backwards look to rho_01 and how we access it -- it's due to Fortran2018's printing of arrays)
   if imag_or_real == 'imag':
      indx = (total_states*total_states) + ( ( jstate - 1 ) * total_states) + istate
   elif imag_or_real == 'real':
      indx = ( ( jstate - 1 ) * total_states) + istate
   
   return indx
   
def animate_bloch_plot( bloch_obj, time_obj, movie_slice_step, filename='bloch_movie.mp4', imagename='bloch_sphere.png' ):
    
   # User options
   zlabelmod = 0.00
   nzlabelmod = -0.12
   timelabelx, timelabely = 0.0, 0.05
   vectorwidth_fact = 2
   pointwidth_fact = 0.125
   time_unit = 'ns'

   # Shortcut for colors
   redcolor = np.array([1.,0.,0.])
   greencolor = np.array([0.,1.,0.])
   bluecolor = np.array([0.,0.,1.])
   yellowcolor = np.array([1.,1.,0.])
   purplecolor = np.array([0.5,0.,0.5])

   # Movie and Path
   b = Bloch()
   b.make_sphere()
   movrange = len(bloch_obj[::movie_slice_step,0])
   xpmov = bloch_obj[::movie_slice_step,0]
   ypmov = bloch_obj[::movie_slice_step,1]
   zpmov = bloch_obj[::movie_slice_step,2]
   tpmov = time_obj[::movie_slice_step]
   check = bloch_obj[::movie_slice_step,0:3]
    
   print("Original number of time points: ",len(time_obj[:]))
   print("Requested slice step: ",movie_slice_step)
   print("Total number of points for movie: ",movrange)
    
   # Bloch sphere labels, vectors, and points
   def plot_setup(obj):
      obj.xlabel = [r'$\mathbf{\hat{x}}$', '']
      obj.ylabel = [r'$\mathbf{\hat{y}}$', '']
      obj.zlabel = [r'$\mathbf{\left|0\right>}$', r'$\mathbf{\left|1\right>}$']
      obj.zlpos = [1.2+zlabelmod,-1.2+nzlabelmod]
      #obj.vector_width = int(vectorwidth_fact*4.)
      obj.point_size = [int(pointwidth_fact*55.)]

   # Save path
   plot_setup(b)
   for i in range(movrange):
      b.add_points(check[i], 'm')

   print("Saving trajectory...")
   b.save(name=imagename)
   print("Trajectory saved!")

   # Movie
   fig = plt.figure()
   ax = fig.add_subplot(111, projection="3d", elev=30, azim=-40)
   sphere = qutip.Bloch(axes=ax)

   def init():
      plot_setup(sphere)
      return ax

   def animate(i):
      sphere.clear()
      sphere.add_vectors([xpmov[i],ypmov[i],zpmov[i]])
      #sphere.add_vectors([0.,0.,zpmov[i]])
      #sphere.add_points([xpmov[:i+1], ypmov[:i+1], zpmov[:i+1]])
      sphere.make_sphere()
      ax.text2D(timelabelx, timelabely,f"{tpmov[i]:0.2f} "+time_unit,transform=ax.transAxes)
      return ax

   ani = animation.FuncAnimation(fig, animate, np.arange(len(xpmov)), init_func=init, blit=False, repeat=False)
                               
   ani.save(filename, fps=60, progress_callback = lambda i, n: printProgressBar(i + 1, n, prefix = 'Movie Progress:', suffix = 'Complete', length = 50))
   
def project_bloch_sphere( rho_0_state_num, rho_1_state_num, rho_total_printed_states, rotframe, f_d, eps1, eps2, rho_filename, print_freqs ):
   # Get column indices for each state
   indx_Re_00 = get_py_indx_rhofile( rho_0_state_num, rho_0_state_num, 'real', rho_total_printed_states )
   indx_Re_01 = get_py_indx_rhofile( rho_0_state_num, rho_1_state_num, 'real', rho_total_printed_states )
   indx_Re_10 = get_py_indx_rhofile( rho_1_state_num, rho_0_state_num, 'real', rho_total_printed_states )
   indx_Re_11 = get_py_indx_rhofile( rho_1_state_num, rho_1_state_num, 'real', rho_total_printed_states )
   indx_Im_00 = get_py_indx_rhofile( rho_0_state_num, rho_0_state_num, 'imag', rho_total_printed_states )
   indx_Im_01 = get_py_indx_rhofile( rho_0_state_num, rho_1_state_num, 'imag', rho_total_printed_states )
   indx_Im_10 = get_py_indx_rhofile( rho_1_state_num, rho_0_state_num, 'imag', rho_total_printed_states )
   indx_Im_11 = get_py_indx_rhofile( rho_1_state_num, rho_1_state_num, 'imag', rho_total_printed_states )

   # Load data
   data = np.loadtxt(rho_filename)
   rho_sub = np.zeros([len(data[:,0]),2,2],dtype=np.cdouble)
   data_new = np.zeros([len(data[:,0]),4])

   # Inform user
   print("Building Bloch sphere where |0> corresponds to eigenstate "+str(rho_0_state_num)+" and |1> corresponds to eigenstate "+str(rho_1_state_num))
   print("For reference, rho_00, rho_01, and rho_10 elements are taken from 'RHO.dat'")
   print("The corresponding column numbers in 'RHO.dat' for the real entries of those elements are: "+str(indx_Re_00)+", "+str(indx_Re_01)+", and "+str(indx_Re_11))

   # Form sub density matrix block
   
   # note for flipped order rho states, i.e. one chooses rho_0_state_num > rho_1_state_num: rho_new = U rho U^dagger, where U = [[0,1],[1,0]]
   # assigning rho_(rho_0_state_num,rho_0_state_num) = b, rho_(rho_1_state_num,rho_1_state_num) = a
   # and rho_(rho_1_state_num,rho_0_state_num) = c,
   # rho originally was [[a,c],[c*, b]], and upon transformation, rho_new = [[b,c*],[c,a]]
   # Notice the flip of c* and c.
   
   rho_sub[:,0,0] = data[:,indx_Re_00] + (1.0j*data[:,indx_Im_00])
   if rho_0_state_num > rho_1_state_num:
      rho_sub[:,0,1] = data[:,indx_Re_01] + (-1.0j*data[:,indx_Im_01])
   else:
      rho_sub[:,0,1] = data[:,indx_Re_01] + (1.0j*data[:,indx_Im_01])
   rho_sub[:,1,0] = np.conjugate(rho_sub[:,0,1])
   rho_sub[:,1,1] = data[:,indx_Re_11] + (1.0j*data[:,indx_Im_11])

   # time
   data_new[:,0] = data[:,0]

   # Change to rotating frame
   if rotframe:
   
      if print_freqs:
         # Get frequencies
         
         # x proj in Bloch sphere = 2 * Re rho_01
         data_new[:,1] = 2.0*np.real(rho_sub[:,0,1])

         # y proj in Bloch sphere = 2 * Im rho_01
         data_new[:,2] = 2.0*np.imag(rho_sub[:,0,1])

         # z proj in Bloch sphere = rho_00 - rho_11
         data_new[:,3] = np.real(rho_sub[:,0,0] - rho_sub[:,1,1])
         
         x_freq = extract_peak_frequency(data_new[:,1], data_new[0,0], data_new[-1,0], len(data_new[:,0]))
         print("Peak frequency of oscillations in the x-axis of the Bloch sphere before rotation transformation is",x_freq,"GHz")
         y_freq = extract_peak_frequency(data_new[:,2], data_new[0,0], data_new[-1,0], len(data_new[:,0]))
         print("Peak frequency of oscillations in the y-axis of the Bloch sphere before rotation transformation is",y_freq,"GHz")
         z_freq = extract_peak_frequency(data_new[:,3], data_new[0,0], data_new[-1,0], len(data_new[:,0]))
         print("Peak frequency of oscillations in the z-axis of the Bloch sphere before rotation transformation is",z_freq,"GHz")
   
      urot = np.zeros((len(data[:,0]),2,2),dtype=np.cdouble)
      urot[:,0,0] = np.exp( -1.0j * 2.0*np.pi * f_d * eps1 / (eps2 - eps1) * data_new[:,0] )
      urot[:,1,1] = np.exp( -1.0j * 2.0*np.pi * f_d * eps2 / (eps2 - eps1) * data_new[:,0] )
      for i in range(len(data[:,0])):
          rho_sub[i,:,:] = np.matmul(urot[i,:,:],rho_sub[i,:,:])
          rho_sub[i,:,:] = np.matmul(rho_sub[i,:,:],np.transpose(np.conjugate(urot[i,:,:])))

   # x proj in Bloch sphere = 2 * Re rho_01
   data_new[:,1] = 2.0*np.real(rho_sub[:,0,1])

   # y proj in Bloch sphere = 2 * Im rho_01
   data_new[:,2] = 2.0*np.imag(rho_sub[:,0,1])

   # z proj in Bloch sphere = rho_00 - rho_11
   data_new[:,3] = np.real(rho_sub[:,0,0] - rho_sub[:,1,1])
   
   if print_freqs:
      if rotframe:
         # Get frequencies
         x_freq = extract_peak_frequency(data_new[:,1], data_new[0,0], data_new[-1,0], len(data_new[:,0]))
         print("Peak frequency of oscillations in the x-axis of the Bloch sphere after rotation transformation is",x_freq,"GHz")
         y_freq = extract_peak_frequency(data_new[:,2], data_new[0,0], data_new[-1,0], len(data_new[:,0]))
         print("Peak frequency of oscillations in the y-axis of the Bloch sphere after rotation transformation is",y_freq,"GHz")
         z_freq = extract_peak_frequency(data_new[:,3], data_new[0,0], data_new[-1,0], len(data_new[:,0]))
         print("Peak frequency of oscillations in the z-axis of the Bloch sphere after rotation transformation is",z_freq,"GHz")
      else:
         # Get frequencies
         x_freq = extract_peak_frequency(data_new[:,1], data_new[0,0], data_new[-1,0], len(data_new[:,0]))
         print("Peak frequency of oscillations in the x-axis of the Bloch sphere is",x_freq,"GHz")
         y_freq = extract_peak_frequency(data_new[:,2], data_new[0,0], data_new[-1,0], len(data_new[:,0]))
         print("Peak frequency of oscillations in the y-axis of the Bloch sphere is",y_freq,"GHz")
         z_freq = extract_peak_frequency(data_new[:,3], data_new[0,0], data_new[-1,0], len(data_new[:,0]))
         print("Peak frequency of oscillations in the z-axis of the Bloch sphere is",z_freq,"GHz")

   # Save new file
   return data_new

   
