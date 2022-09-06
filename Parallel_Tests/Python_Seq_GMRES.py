#Named Sequential, but utilises all cores on CPU
from scipy.io import mmread
from scipy.sparse.linalg import *
import numpy as np
import time
A = mmread('thermal1.mtx')
b = mmread('thermal1_b.mtx')
tols = [1e-02, 1e-03, 1e-04, 1e-05]
for tol in tols:
    print("----tol="+str(tol)+"----")
    start = time.time()
    x, exitCode = gmres(A, b, tol=tol, atol=0, restart = 50)
    end = time.time()
    print("Time s: "+str((end-start)))
    print("|Ax-b| = "+str(np.linalg.norm(A.dot(x)-b[:,0])))