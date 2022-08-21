#Analysis
from Solvers import *
from Examples import *
from firedrake import *

def Example_1_Analysis(al = 0, m = 1):
    a = Example_1(al = al, m=m)
    u_e = a[5]
    u_h = Solve_Limit_Method(*a[0:5], Order = 3)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("LM eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("LM eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_MMAP_Method(*a[0:5], Order = 3)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_Singular_Pert_Method(*a[0:5], Order = 3)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("SPM eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("SPM eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_AP_Method(*a[0:5], Order = 3)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("AP eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("AP eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_PF_Method(*a[0:5], Order = 3)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("PF eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("PF eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)
    
    '''
    #DN Section
    u_h = Solve_DN_Method(*a[0:5], Order = 3, eps_0 = 1e-3)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("DN eps 1e-15 0 1e-3 ||u - u_exact||_L2: %.1e" % L2_error)
    print("DN eps 1e-15 0 1e-3 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_DN_Method(*a[0:5], Order = 3, eps_0 = 1e-6)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("DN eps 1e-15 0 1e-6 ||u - u_exact||_L2: %.1e" % L2_error)
    print("DN eps 1e-15 0 1e-6 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_DN_Method(*a[0:5], Order = 3, eps_0 = 1e-9)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("DN eps 1e-15 0 1e-9 ||u - u_exact||_L2: %.1e" % L2_error)
    print("DN eps 1e-15 0 1e-9 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_DN_Method(*a[0:5], Order = 3, eps_0 = 1e-12)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("DN eps 1e-15 0 1e-12 ||u - u_exact||_L2: %.1e" % L2_error)
    print("DN eps 1e-15 0 1e-12 ||u - u_exact||_H1: %.1e" % H1_error)
    '''
    
    a = Example_1(eps = 0.1, al = al, m=m)
    u_e = a[5]
    u_h = Solve_Limit_Method(*a[0:5], eps = 0.1, Order = 3)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("LM eps 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("LM eps 0.1 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_MMAP_Method(*a[0:5], eps = 0.1, Order = 3)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP eps 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP eps 0.1 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_Singular_Pert_Method(*a[0:5], eps = 0.1, Order = 3)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("SPM eps 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("SPM eps 0.1 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_AP_Method(*a[0:5], eps = 0.1, Order = 3)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("AP eps 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("AP eps 0.1 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_PF_Method(*a[0:5], eps = 0.1, Order = 3)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("PF eps 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("PF eps 0.1 ||u - u_exact||_H1: %.1e" % H1_error)

def Example_2_Mesh_Comp(Order = 3):
    E2 = Example_2()
    u_e = E2[5]
    eps = 1e-15
    u_h = Solve_Limit_Method(*E2, Order = Order, eps = eps)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("LM ||u - u_exact||_L2: %.1e" % L2_error)
    print("LM ||u - u_exact||_H1: %.1e" % H1_error) 
    
    u_h = Solve_Limit_STAB_Method(*E2, Order = Order, eps = eps)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("LM ||u - u_exact||_L2: %.1e" % L2_error)
    print("LM ||u - u_exact||_H1: %.1e" % H1_error) 
    
    E2 = Example_2(Restrict_CL = True)
    u_e = E2[5]
    eps = 1e-15
    u_h = Solve_Limit_Method(*E2, Order = Order, eps = eps)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("LM CL ||u - u_exact||_L2: %.1e" % L2_error)
    print("LM CL ||u - u_exact||_H1: %.1e" % H1_error)    
    
def Example_2_Analysis(a = 0.05, Order = 3, Size = 40):
    E2 = Example_2(a = a, Size = Size)
    u_e = E2[5]
    eps = 1e-15
    u_h = Solve_Limit_Method(*E2, Order = Order, eps = eps)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("LM eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("LM eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_MMAP_Method(*E2, Order = Order, eps = eps)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_MMAP_STAB_Method(*E2, Order = Order, eps = eps)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP STAB eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP STAB eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_Singular_Pert_Method(*E2, Order = Order, eps = eps)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("SPM eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("SPM eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_AP_Method(*E2, Order = Order, eps = eps)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("AP eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("AP eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_PF_Method(*E2, Order = Order, eps = eps)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("PF eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("PF eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)
    
    
    
    '''
    #DN Section
    u_h = Solve_DN_Method(*a[0:5], Order = 3, eps_0 = 1e-3)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("DN eps 1e-15 0 1e-3 ||u - u_exact||_L2: %.1e" % L2_error)
    print("DN eps 1e-15 0 1e-3 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_DN_Method(*a[0:5], Order = 3, eps_0 = 1e-6)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("DN eps 1e-15 0 1e-6 ||u - u_exact||_L2: %.1e" % L2_error)
    print("DN eps 1e-15 0 1e-6 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_DN_Method(*a[0:5], Order = 3, eps_0 = 1e-9)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("DN eps 1e-15 0 1e-9 ||u - u_exact||_L2: %.1e" % L2_error)
    print("DN eps 1e-15 0 1e-9 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_DN_Method(*a[0:5], Order = 3, eps_0 = 1e-12)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("DN eps 1e-15 0 1e-12 ||u - u_exact||_L2: %.1e" % L2_error)
    print("DN eps 1e-15 0 1e-12 ||u - u_exact||_H1: %.1e" % H1_error)
    '''
    
    E2 = Example_2(eps = 0.1, a = a, Size = Size)
    u_e = E2[5]
    u_h = Solve_Limit_Method(*E2, eps = 0.1, Order = Order)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("LM eps 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("LM eps 0.1 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_MMAP_Method(*E2, eps = 0.1, Order = Order)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP eps 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP eps 0.1 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_MMAP_STAB_Method(*E2, eps = 0.1, Order = Order)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP STAB eps 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP STAB eps 0.1 ||u - u_exact||_H1: %.1e" % H1_error)   
    
    u_h = Solve_Singular_Pert_Method(*E2, eps = 0.1, Order = Order)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("SPM eps 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("SPM eps 0.1 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_AP_Method(*E2, eps = 0.1, Order = Order)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("AP eps 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("AP eps 0.1 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_PF_Method(*E2, eps = 0.1, Order = Order)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("PF eps 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("PF eps 0.1 ||u - u_exact||_H1: %.1e" % H1_error)
    
def Example_2_Analysis_Stab(a = 0.05, Order = 3, Size = 40):
    E2 = Example_2(a = a, Size = Size)
    u_e = E2[5]
    eps = 1e-15

    u_h = Solve_MMAP_STAB_Method(*E2, Order = Order, eps = eps, sigma = 1.)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP STAB eps 1e-15 sigma 1. ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP STAB eps 1e-15 sigma 1. ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_MMAP_STAB_Method(*E2, Order = Order, eps = eps, sigma = 0.1)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP STAB eps 1e-15 sigma 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP STAB eps 1e-15 sigma 0.1 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_MMAP_STAB_Method(*E2, Order = Order, eps = eps, sigma = 0.01)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP STAB eps 1e-15 sigma 0.01 ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP STAB eps 1e-15 sigma 0.01 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_MMAP_STAB_Method(*E2, Order = Order, eps = eps, sigma = 0.001)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP STAB eps 1e-15 sigma 0.001 ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP STAB eps 1e-15 sigma 0.001 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_MMAP_STAB_Method(*E2, Order = Order, eps = eps, sigma = 0.0001)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP STAB eps 1e-15 sigma 0.0001 ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP STAB eps 1e-15 sigma 0.0001 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_MMAP_Method(*E2, Order = Order, eps = eps)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)
    
def Example_2_Analysis_PF_Stab(a = 0.05, Order = 3, Size = 40):
    E2 = Example_2(a = a, Size = Size)
    u_e = E2[5]
    eps = 1e-15

    u_h = Solve_PF_STAB_Method(*E2, Order = Order, eps = eps, sigma = 1.)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP STAB eps 1e-15 sigma 1. ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP STAB eps 1e-15 sigma 1. ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_PF_STAB_Method(*E2, Order = Order, eps = eps, sigma = 0.1)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP STAB eps 1e-15 sigma 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP STAB eps 1e-15 sigma 0.1 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_PF_STAB_Method(*E2, Order = Order, eps = eps, sigma = 0.01)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP STAB eps 1e-15 sigma 0.01 ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP STAB eps 1e-15 sigma 0.01 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_PF_STAB_Method(*E2, Order = Order, eps = eps, sigma = 0.001)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP STAB eps 1e-15 sigma 0.001 ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP STAB eps 1e-15 sigma 0.001 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_PF_STAB_Method(*E2, Order = Order, eps = eps, sigma = 0.0001)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP STAB eps 1e-15 sigma 0.0001 ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP STAB eps 1e-15 sigma 0.0001 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_PF_Method(*E2, Order = Order, eps = eps)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)
 
    
def Example_3_Analysis():
    a = Example_3(eps = 1e-7)
    u_e = a[5]
    
    u_h = Solve_Limit_Method(*a, Order = 2, eps = 1e-7)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("LM eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("LM eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_PF_Method(*a, Order = 2, eps = 1e-7)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("PF eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("PF eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_PF_STAB_Method(*a, Order = 2, eps = 1e-7)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("PF_STAB eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("PF_STAB eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)    

        
    a = Example_3(eps = 0.1)
    u_e = a[5]
    
    u_h = Solve_PF_Method(*a, eps = 0.1, Order = 2)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("PF eps 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("PF eps 0.1 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_PF_STAB_Method(*a, eps = 0.1, Order = 2)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("PF_STAB eps 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("PF_STAB eps 0.1 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_MMAP_Method(*a, eps = 0.1, Order = 2)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP eps 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP eps 0.1 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_Singular_Pert_Method(*a, eps = 0.1, Order = 2)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("SPM eps 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("SPM eps 0.1 ||u - u_exact||_H1: %.1e" % H1_error)
    
def Example_4_Analysis():
    a = Example_4(eps = 1e-7)
    u_e = a[5]
    
    '''
    u_h = Solve_PF_Method(*a, Order = 2, eps = 1e-7)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("PF eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("PF eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)
    '''
    u_h = Solve_PF_STAB_Method(*a, Order = 2, eps = 1e-7)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("PF eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("PF eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)    

        
    a = Example_4(eps = 0.1)
    u_e = a[5]
    
    u_h = Solve_PF_Method(*a, eps = 0.1, Order = 2)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("PF eps 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("PF eps 0.1 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_PF_STAB_Method(*a, eps = 0.1, Order = 2)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("PF eps 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("PF eps 0.1 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_MMAP_Method(*a, eps = 0.1, Order = 2)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("MMAP eps 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("MMAP eps 0.1 ||u - u_exact||_H1: %.1e" % H1_error)
    
    u_h = Solve_Singular_Pert_Method(*a, eps = 0.1, Order = 2)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("SPM eps 0.1 ||u - u_exact||_L2: %.1e" % L2_error)
    print("SPM eps 0.1 ||u - u_exact||_H1: %.1e" % H1_error)
    
    