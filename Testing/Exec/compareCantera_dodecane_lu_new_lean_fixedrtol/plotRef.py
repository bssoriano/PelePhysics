import numpy as np
import sys
sys.path.append('utils')
import fileio as io
import matplotlib.pyplot as plt
from plotsUtil import *
import argparse


Aref = np.load('ref/Aref.npy')
fig = plt.figure()
plt.plot(Aref[:,0], Aref[:,1], color='k', linewidth=3)
prettyLabels("Time [s]", "T [K]", 14)
plt.savefig("reference_sol.png")
plt.savefig("reference_sol.eps")
