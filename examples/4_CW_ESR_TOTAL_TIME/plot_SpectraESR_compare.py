from helper import *
from os import sys

print(sys.argv)
numpoints = int(sys.argv[1])
filename = sys.argv[1+numpoints+1]
outputname = sys.argv[1+numpoints+2]
datafiles, names = [], []
for i in range(numpoints):
   foldername = sys.argv[2+i]
   datafiles.append(foldername + '/' + filename)
   names.append(foldername)
plot_esr_compare( freq_diff=0.1, datafiles=datafiles, names=names, savename=outputname )
