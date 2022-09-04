from Solvers import *
from Examples import *
from firedrake import *

E1 = Example_1(al = 2, m=10, eps = 1e-10, Size = 160)
u_e = E1[5]
u_h = Solve_MMAP_Method(*E1, Order = 2, eps = 1e-10, save = True)
#Save creates a PVD file with u_h, f, ... in output_graph/
L2_error = norm(u_h - u_e, "L2")
H1_error = norm(u_h - u_e, "H1")
print("MMAP eps 1e-10 ||u - u_exact||_L2: %.1e" % L2_error)
print("MMAP eps 1e-10 ||u - u_exact||_H1: %.1e" % H1_error)