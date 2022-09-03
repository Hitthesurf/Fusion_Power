#makePics
from Solvers import *
from Examples import *
from firedrake import *
from decimal import Decimal
import matplotlib.pyplot as plt
import numpy as np


def Map_f():
    E4 = Example_4(eps = 1e-15, Size = 40)
    epss = 10**(np.arange(-15,1,1,dtype=float))
    u_e = E4[5]
    for eps in epss:
        u_h = Solve_MMAP_STAB_Method(*E4, Order = 2, eps = eps)
        print("eps: "+str(eps) +", L2: %.5f" % norm(u_e-u_h, "L2"))
        print("eps: "+str(eps) +", H1: %.5f" % norm(u_e-u_h, "H1"))

def GenLH(Example, Solver, epss):
    L2 = []
    H1 = []
    for eps in epss:
        a = Example(eps = eps)
        u_e = a[5]
        u_h = Solver(*a, Order = 2, eps = eps)
        L2.append(norm(u_h - u_e, "L2"))
        H1.append(norm(u_h - u_e, "H1"))
    return L2, H1
    
def GenSize(Example, Solver, Sizes):
    L2 = []
    H1 = []
    for Size in Sizes:
        a = Example(Size = Size)
        u_e = a[5]
        u_h = Solver(*a, Order = 2, eps = 1e-10)
        L2.append(norm(u_h - u_e, "L2"))
        H1.append(norm(u_h - u_e, "H1"))
    return L2, H1
    
def PlotFig(Example, Methods, Names, epss, strID):
    #Create Figure
    H1_Over = []
    plt.figure(figsize=(10, 10), dpi=200)     
    for ID in range(len(Methods)):
        Name = Names[ID]
        Method = Methods[ID]
        L2, H1 = GenLH(Example, Method, epss)
        plt.plot(epss, L2, label = Name)
        plt.scatter(epss, L2, color = 'red', zorder = 2)
        H1_Over.append(H1)   
    plt.xscale('log')
    plt.yscale('log')
    #Save Figure
    plt.xlabel(r'$\varepsilon$')
    plt.ylabel(r'L2 Error')
    plt.rc('font', size=16)
    plt.legend()
    plt.savefig("Pics/LHSims/"+strID+"L2.png")
    
    #Create for H1 figure
    plt.figure(figsize=(10, 10), dpi=200) 
    
    for ID in range(len(Methods)):
        Name = Names[ID]
        Method = Methods[ID]
        plt.plot(epss, H1_Over[ID], label = Name) 
        plt.scatter(epss, H1_Over[ID], color = 'red', zorder = 2)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$\varepsilon$')
    plt.ylabel(r'H1 Error')
    plt.rc('font', size=16)
    #Save Figure
    plt.legend()
    plt.savefig("Pics/LHSims/"+strID+"H1.png")
    
def Plot1():
    # Size = 20
    #Over head same for every example
    Methods = [Solve_MMAP_Method,
               Solve_AP_Method, Solve_PF_Method] #, Solve_Limit_Method
    Method_Names = ["MMAP", "AP", "PF"]
    epss = 10**(np.arange(-15,1,1,dtype=float))
    
    
    E1a = lambda eps: Example_1(eps = eps, Size = 20)
    PlotFig(E1a, Methods, Method_Names, epss, "E1a_MMAP_AP_PF")
    E1b = lambda eps: Example_1(al = 2, m = 1, eps = eps, Size = 20)
    PlotFig(E1b, Methods, Method_Names, epss, "E1b_MMAP_AP_PF")
    E1c = lambda eps: Example_1(al = 2, m = 10, eps = eps, Size = 20)
    PlotFig(E1c, Methods, Method_Names, epss, "E1c_MMAP_AP_PF")
    
    #Methods = [Solve_Limit_Method, Solve_Singular_Pert_Method]

def GenTable(Example, Methods, Names, Sizes, strID):
    for ID in range(len(Methods)):
        Method = Methods[ID]
        Name = Names[ID]
        L2, H1 = GenSize(Example, Method, Sizes)
        print("-----------StartInfo----------")
        print("Method: " + Name + ", strID: " + strID)
        print("Sizes: "+str(Sizes))
        print("L2: " + str(L2))
        print("H1: " + str(H1))
        print("-----------ENDInfo------------")
        

def Plot2():
    # Size = 20
    #Over head same for every example
    Methods = [Solve_MMAP_Method,
               Solve_Limit_Method, Solve_Singular_Pert_Method] #, Solve_Limit_Method
    Method_Names = ["MMAP", "LM", "SP"]
    epss = 10**(np.arange(-15,1,1,dtype=float))
    
    
    E1a = lambda eps: Example_1(eps = eps, Size = 20)
    PlotFig(E1a, Methods, Method_Names, epss, "E1a_MMAP_LM_SP")
    E1b = lambda eps: Example_1(al = 2, m = 1, eps = eps, Size = 20)
    PlotFig(E1b, Methods, Method_Names, epss, "E1b_MMAP_LM_SP")
    E1c = lambda eps: Example_1(al = 2, m = 10, eps = eps, Size = 20)
    PlotFig(E1c, Methods, Method_Names, epss, "E1c_MMAP_LM_SP") 

def PlotE2():

    epss = 10**(np.arange(-15,1,1,dtype=float))
    # Size Fixxed by mesh
    Methods = [Solve_MMAP_Method, Solve_PF_Method]
    Method_Names = ["MMAP", "PF"]
    
       
    E2 = lambda eps: Example_2(eps = eps, Restrict_CL = False)
    PlotFig(E2, Methods, Method_Names, epss, "E2/E2_Normal")
    
    #Gamma_In
    E2R = lambda eps: Example_2(eps = eps, Restrict_CL = True)
    PlotFig(E2R, Methods, Method_Names, epss, "E2/E2_IN")
    
    #STAB, sigma = 0.1
    Methods = [Solve_MMAP_STAB_Method, Solve_PF_STAB_Method]
    Method_Names = ["MMAP_STAB", "PF_STAB"]
    PlotFig(E2, Methods, Method_Names, epss, "E2/E2_STAB")
    
def PlotE4():
    epss = 10**(np.arange(-15,1,1,dtype=float))
    #Size 20 x 20
    Methods = [Solve_MMAP_STAB_Method, Solve_PF_STAB_Method]
    Method_Names = ["MMAP_STAB", "PF_STAB"]
    
    E4 = lambda eps: Example_4(eps = eps, Size = 20)
    PlotFig(E4, Methods, Method_Names, epss, "E4/E4_STAB")
    
def PlotE5():
    #Size fixxed by mesh
    epss = 10**(np.arange(-10,1,1,dtype=float))
    Methods = [Solve_PF_STAB_Method, Solve_MMAP_STAB_Method]
    Method_Names = ["PF_STAB", "MMAP_STAB"]
    E5 = lambda eps: Example_5(eps=eps, Restrict_CL = False)
    PlotFig(E5, Methods, Method_Names, epss, "E5/E5_STAB")
    

def GetPVD_Diff_E3():
    #Size 80x80
    #dof 51842
    eps = 1e-10
    
    #PF_Stab
    E3 = Example_3(Size = 80, eps = eps)
    u_e = E3[5]
    
    u_h = Solve_PF_STAB_Method(*E3, eps = eps, Order = 2, save = True)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("SPM eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("SPM eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)
    
    #MMAP_Stab
    u_h = Solve_MMAP_STAB_Method(*E3, eps = eps, Order = 2, save = True)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("SPM eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("SPM eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)
    
def GetPVD_Diff_E4():
    #Size 80x80
    #dof 51842
    eps = 1e-10
    
    #PF_Stab
    E3 = Example_4(Size = 80, eps = eps)
    u_e = E3[5]
    u_h = Solve_PF_STAB_Method(*E3, eps = eps, Order = 2, save = True)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("SPM eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("SPM eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)
    
    #MMAP_Stab
    u_h = Solve_MMAP_STAB_Method(*E3, eps = eps, Order = 2, save = True)
    L2_error = norm(u_h - u_e, "L2")
    H1_error = norm(u_h - u_e, "H1")
    print("SPM eps 1e-15 ||u - u_exact||_L2: %.1e" % L2_error)
    print("SPM eps 1e-15 ||u - u_exact||_H1: %.1e" % H1_error)   
    
    

def Table1():
    Sizes = [10,20,40, 80, 160]
    eps = 1e-10
    Methods = [Solve_PF_Method, Solve_MMAP_Method]
    Method_Names = ["PF", "MMAP"]
    
    E1a = lambda Size: Example_1(eps = eps, Size = Size)
    GenTable(E1a, Methods, Method_Names, Sizes, "E1a_MMAP_PF")
    
    E1b = lambda Size: Example_1(al=2, m=1, eps=eps, Size = Size)
    GenTable(E1b, Methods, Method_Names, Sizes, "E1b_MMAP_PF")
    
    E1c = lambda Size: Example_1(al=2, m=10, eps=eps, Size = Size)
    GenTable(E1c, Methods, Method_Names, Sizes, "E1c_MMAP_PF")
    
def TableE4():
    Sizes = [10,20,40, 80, 160]
    eps = 1e-10
    #a=0.05
    Methods = [Solve_MMAP_STAB_Method, Solve_PF_STAB_Method]
    Method_Names = ["MMAP_STAB", "PF_STAB"]
   
    E4 = lambda Size: Example_4(eps = eps, Size = Size)
    GenTable(E4, Methods, Method_Names, Sizes, "E4_STAB")

def UpdateSave(to_save, a, ID):
    f = a[0]
    mesh = a[2]
    V = FunctionSpace(mesh, 'CG', 4)
    u = Function(V).interpolate(a[5])
    f = Function(V).interpolate(a[0])
    u.rename('u_' + ID)
    f.rename('f_' + ID)
    to_save.append(u)
    to_save.append(f)
        #Save file
    File("PVD/FU_"+ID+".pvd").write(u, f) 
    return to_save

def Getfu():
    #Create pvd file
    #E1a_eps_1
    to_save = []
    ID = 'E1a_eps_1'
    a = Example_1(al = 0, m=1, eps = 0.1)
    to_save = UpdateSave(to_save, a, ID)
    
    #E1a_eps_10
    ID = 'E1a_eps_10'
    a = Example_1(al = 0, m=1, eps = 1e-10)
    to_save = UpdateSave(to_save, a, ID)
    
    #E1b_eps_1
    ID = 'E1b_eps_1'
    a = Example_1(al = 2, m=1, eps = 0.1)
    to_save = UpdateSave(to_save, a, ID)
    
    #E1b_eps_10
    ID = 'E1b_eps_10'
    a = Example_1(al = 2, m=1, eps = 1e-10)
    to_save = UpdateSave(to_save, a, ID)
    
    #E1c_eps_1
    ID = 'E1c_eps_1'
    a = Example_1(al = 2, m=10, eps = 0.1)
    to_save = UpdateSave(to_save, a, ID)
    
    #E1c_eps_10
    ID = 'E1c_eps_10'
    a = Example_1(al = 2, m=10, eps = 1e-10)
    to_save = UpdateSave(to_save, a, ID)
    
    #E2_eps_1
    ID = 'E2_eps_1'
    a = Example_2(eps = 0.1)
    to_save = UpdateSave(to_save, a, ID)
    
    #E2_eps_10
    ID = 'E2_eps_10'
    a = Example_2(eps = 1e-10)
    to_save = UpdateSave(to_save, a, ID) 
    
    #E3_eps_1
    ID = 'E3_eps_1'
    a = Example_3(eps = 0.1)
    to_save = UpdateSave(to_save, a, ID)
    
    #E3_eps_10
    ID = 'E3_eps_10'
    a = Example_3(eps = 1e-10)
    to_save = UpdateSave(to_save, a, ID) 
    
    #E4_eps_1
    ID = 'E4_eps_1'
    a = Example_4(eps = 0.1)
    to_save = UpdateSave(to_save, a, ID)
    
    #E4_eps_10
    ID = 'E4_eps_10'
    a = Example_4(eps = 1e-10)
    to_save = UpdateSave(to_save, a, ID) 
    
    #    #E5_eps_1
    ID = 'E5_eps_1'
    a = Example_5(eps = 0.1)
    to_save = UpdateSave(to_save, a, ID)
    
    #E5_eps_10
    ID = 'E5_eps_10'
    a = Example_5(eps = 1e-10)
    to_save = UpdateSave(to_save, a, ID) 
    

    