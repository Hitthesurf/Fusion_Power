E4 = Example_4(eps = 1e-15, Size = 40)
epss = 10**(np.arange(-15,1,1,dtype=float))
u_e = E4[5]
for eps in epss:
    u_h = Solve_MMAP_STAB_Method(*E4, Order = 2, eps = eps)
    print("eps: "+str(eps) +", L2: %.5f" % norm(u_e-u_h, "L2"))
    print("eps: "+str(eps) +", H1: %.5f" % norm(u_e-u_h, "H1"))