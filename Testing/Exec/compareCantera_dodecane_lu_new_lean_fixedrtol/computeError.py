import numpy as np
import sys
sys.path.append('utils')
import fileio as io
import matplotlib.pyplot as plt
from plotsUtil import *
import argparse

def computeFieldErr(ind, A, Aref):
    interpField = np.interp(A[:,0], Aref[:,0], Aref[:,ind])
    if A[-1,0] < 0.099:
        error = np.nan
    else:
        error = np.linalg.norm(A[:,ind] - interpField) / A.shape[0]
    return error

def computeError(file, fileRef_np=None, fileRef= 'ref/PPreaction.txt'):
    if fileRef_np is not None:
        Aref = np.load(fileRef_np)
    else:
        Aref = io.readMultiColFile(fileRef)
    A = io.readMultiColFile(file)
    return computeFieldErr(1, A, Aref)


parser = argparse.ArgumentParser(description="Error comp")
parser.add_argument(
    "-f",
    "--result_file",
    type=str,
    metavar="",
    required=True,
    help="File with results of 0D",
    default=None,
)
parser.add_argument(
    "-l",
    "--log_file",
    type=str,
    metavar="",
    required=True,
    help="accuracy log",
    default=None,
)
args = parser.parse_args()

err=computeError(args.result_file, fileRef_np='ref/Aref.npy', fileRef= 'ref/PPreaction.txt')
f = open(args.log_file, 'w+')
f.write(f"{err}\n")
f.close()
