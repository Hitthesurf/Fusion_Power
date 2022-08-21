#Solvers
from firedrake import *
import numpy as np

def SAVE(to_save, f, b, mesh, u_h, u_e, Order, name, Extra):
    path_name_file = "output_graph/"+name+"_solution.pvd"
    U = FunctionSpace(mesh, "CG", Order)
    vecU = VectorFunctionSpace(mesh, "CG", Order)
    b_h = Function(vecU).interpolate(b)
    f_h = Function(U).interpolate(f)
    f_h.rename("f_h")
    b_h.rename("b_h")  
    to_save.append(f_h)
    to_save.append(b_h)
    if Extra:
        u_e_h = Function(U).interpolate(u_e)
        the_error = Function(U).interpolate(u_h-u_e_h)
        grad_u_e_h = Function(vecU).interpolate(grad(u_h))
        u_e_h.rename("u_e_h")
        the_error.rename("the_error")
        grad_u_e_h.rename("grad_u_e_h")
        to_save.append(u_e_h)
        to_save.append(the_error)
        to_save.append(grad_u_e_h)
    File(path_name_file).write(*to_save) 

def Solve_MMAP_Method(f, b, mesh, dOmegaD, dOmegaIN, dOmegaDVal=0,
                      A_para = 1, A_perp = None, Order = 1, eps = 1e-15,
                      save = False, BC_As_Exact = True):
    V = FunctionSpace(mesh, "CG", Order)
    L = FunctionSpace(mesh, "CG", Order)   
    Z = V*L
    
    print(f"MMAP, Z.dim() = {Z.dim()}, Order = {Order}, eps = {eps}.")
    if (eps==0):
        print("LM, As eps=0.")
        
    
    I = Identity(mesh.geometric_dimension())
    if A_perp == None:
        A_perp = I
    
    #Set boundary conditions
    bc0 = DirichletBC(Z.sub(0), dOmegaDVal, dOmegaD)
    bc1 = DirichletBC(Z.sub(1), dOmegaDVal, dOmegaD)
    bcs = [bc0, bc1]
    if (len(dOmegaIN)>0):
        bc2 = DirichletBC(Z.sub(1), 0, dOmegaIN)
        bcs.append(bc2)
    
    grad_para = lambda alpha: dot(outer(b,b), grad(alpha))
    grad_perp = lambda alpha: dot(I-outer(b,b), grad(alpha))
    a_perp = lambda alpha, beta: inner(A_perp*grad_perp(alpha), grad_perp(beta))
    a_para = lambda alpha, beta: inner(A_para*grad_para(alpha), grad_para(beta))    
    
    z = Function(Z)
    u, q = split(z)
    u_, q_ = split(TestFunction(Z))

    F = (
        + a_perp(u, u_)*dx
        + a_para(q, u_)*dx
        - inner(f, u_)*dx
        + a_para(u, q_)*dx        
        )  
    
    if eps>0:
        F = (
            + a_perp(u, u_)*dx
            + a_para(q, u_)*dx
            - inner(f, u_)*dx
            + a_para(u, q_)*dx
            - eps*a_para(q, q_)*dx
            )   
    
    solve(F==0,z,bcs)    
    (u_h, q_h) = z.split()
    u_h.rename("u_h")
    q_h.rename("q_h") 
    
    #Save Data
    if save:
        to_save = [u_h, q_h]
        name = "MMAP"
        if eps == 0:
            name = "LM"
        SAVE(to_save, f, b, mesh, u_h, dOmegaDVal, Order, name, BC_As_Exact)
    return u_h   

    

def Solve_MMAP_STAB_Method(f, b, mesh, dOmegaD, dOmegaIN, dOmegaDVal=0,
                           A_para = 1, A_perp = None, Order = 1, eps = 1e-15,
                           sigma = 0.1, save = False, BC_As_Exact = True):
    V = FunctionSpace(mesh, "CG", Order) 
    Z = V*V
    
    print(f"MMAP_STAB, Z.dim() = {Z.dim()}, Order = {Order}, eps = {eps}.")
    if (eps==0):
        print("LM_STAB, As eps=0.")
    I = Identity(mesh.geometric_dimension())
    if A_perp == None:
        A_perp = I
    
    #Set boundary conditions
    bc0 = DirichletBC(Z.sub(0), dOmegaDVal, dOmegaD)
    bc1 = DirichletBC(Z.sub(1), dOmegaDVal, dOmegaD)
    bcs = [bc0, bc1]
    
    grad_para = lambda alpha: dot(outer(b,b), grad(alpha))
    grad_perp = lambda alpha: dot(I-outer(b,b), grad(alpha))
    a_perp = lambda alpha, beta: inner(A_perp*grad_perp(alpha), grad_perp(beta))
    a_para = lambda alpha, beta: inner(A_para*grad_para(alpha), grad_para(beta))    
    
    z = Function(Z)
    u, q = split(z)
    u_, q_ = split(TestFunction(Z))
    
    
    
    F = (
        + a_perp(u, u_)*dx
        + a_para(q, u_)*dx
        - inner(f, u_)*dx
        + a_para(u, q_)*dx
        - sigma*inner(q,q_)*dx
        )  
    
    if eps>0:
        F = (
            + a_perp(u, u_)*dx
            + a_para(q, u_)*dx
            - inner(f, u_)*dx
            + a_para(u, q_)*dx
            - eps*a_para(q, q_)*dx
            - sigma*inner(q,q_)*dx
            ) 

    solve(F==0,z,bcs)    
    (u_h, q_h) = z.split()
    u_h.rename("u_h")
    q_h.rename("q_h")    
    if save:
        to_save = [u_h, q_h]
        name = "MMAP_STAB"
        if eps == 0:
            name = "LM_STAB"
        SAVE(to_save, f, b, mesh, u_h, dOmegaDVal, Order, name, BC_As_Exact)
    return u_h    


def Solve_Limit_Method(f, b, mesh, dOmegaD, dOmegaIN, dOmegaDVal=0, 
                       A_para = 1, A_perp = None, Order = 1, eps = None,
                       save = False, BC_As_Exact = True):
    return Solve_MMAP_Method(f, b, mesh, dOmegaD, dOmegaIN, dOmegaDVal=dOmegaDVal, 
                             A_para = A_para, A_perp = A_perp, Order = Order, eps = 0,
                             save = save, BC_As_Exact = BC_As_Exact)

def Solve_Limit_STAB_Method(f, b, mesh, dOmegaD, dOmegaIN, dOmegaDVal=0, 
                       A_para = 1, A_perp = None, Order = 1, eps = None,
                            sigma = 0.1, save = False, BC_As_Exact = True):
    return Solve_MMAP_STAB_Method(f, b, mesh, dOmegaD, dOmegaIN, dOmegaDVal=dOmegaDVal, 
                             A_para = A_para, A_perp = A_perp, Order = Order, eps = 0,
                                  sigma = sigma, save = save, BC_As_Exact = BC_As_Exact)

   

def Solve_Singular_Pert_Method(f, b, mesh, dOmegaD, dOmegaIN, dOmegaDVal=0,
                               A_para = 1, A_perp = None, Order = 1, eps = 1e-15,
                               save = False, BC_As_Exact = True):
    U = FunctionSpace(mesh, "CG", Order)
    
    print(f"SP, Z.dim() = {U.dim()}, Order = {Order}, eps = {eps}.")
    I = Identity(mesh.geometric_dimension())
    if A_perp == None:
        A_perp = I
        
    #Set boundary conditions
    bc0 = DirichletBC(U, dOmegaDVal, dOmegaD)
    bcs = [bc0] 
    
    grad_para = lambda alpha: dot(outer(b,b), grad(alpha))
    grad_perp = lambda alpha: dot(I-outer(b,b), grad(alpha))
    a_perp = lambda alpha, beta: inner(A_perp*grad_perp(alpha), grad_perp(beta))
    a_para = lambda alpha, beta: inner(A_para*grad_para(alpha), grad_para(beta))
    
    u = TrialFunction(U)
    psi = TestFunction(U)
    u_h = Function(U)
    
    a = ((1/eps)*a_para(u,psi) + a_perp(u, psi))*dx
    F = inner(f, psi)*dx
    
    solve(a==F, u_h, bcs)
    u_h.rename("u_h")
        
    #Save Data
    if save:
        to_save = [u_h]
        SAVE(to_save, f, b, mesh, u_h, dOmegaDVal, Order, "SP", BC_As_Exact)  
    return u_h

def Solve_AP_Method(f, b, mesh, dOmegaD, dOmegaIN, dOmegaDVal=0,
                    A_para = 1, A_perp = None, Order = 1, eps = 1e-15,
                    save = False, BC_As_Exact = True):
    V = FunctionSpace(mesh, "CG", Order)
    L = FunctionSpace(mesh, "CG", Order)   
    Z = V*L*V*V*L
    
    print(f"AP, Z.dim() = {Z.dim()}, Order = {Order}, eps = {eps}.")
    I = Identity(mesh.geometric_dimension())
    if A_perp == None:
        A_perp = I
    
    #Set boundary conditions
    bc0 = DirichletBC(Z.sub(0), dOmegaDVal, dOmegaD)
    bc1 = DirichletBC(Z.sub(1), dOmegaDVal, dOmegaD)
    
    bc3 = DirichletBC(Z.sub(2), dOmegaDVal, dOmegaD)
    bc4 = DirichletBC(Z.sub(3), dOmegaDVal, dOmegaD)
    bc5 = DirichletBC(Z.sub(4), dOmegaDVal, dOmegaD)
    
    bcs = [bc0, bc1, bc3, bc4, bc5]
    if (len(dOmegaIN)>0):
        bc2 = DirichletBC(Z.sub(1), 0, dOmegaIN)
        bc6 = DirichletBC(Z.sub(4), 0, dOmegaIN)
        bcs.append(bc2)
        bcs.append(bc6)
    
    grad_para = lambda alpha: dot(outer(b,b), grad(alpha))
    grad_perp = lambda alpha: dot(I-outer(b,b), grad(alpha))
    a_perp = lambda alpha, beta: inner(A_perp*grad_perp(alpha), grad_perp(beta))
    a_para = lambda alpha, beta: inner(A_para*grad_para(alpha), grad_para(beta))    
    
    z = Function(Z)
    p, lam, q, l, mu = split(z)
    eta, kappa, xi, chi, tau = split(TestFunction(Z))
    
    F = (
        + a_perp(p, eta)*dx #1
        + a_perp(q, eta)*dx 
        + a_para(eta, lam)*dx
        - inner(f, eta)*dx
        + a_para(p, kappa)*dx #2
        + a_para(q, xi)*dx #3
        + eps*a_perp(q, xi)*dx
        + eps*a_perp(p, xi)*dx
        + inner(l, xi)*dx
        - eps*inner(f, xi)*dx
        + inner(q, chi)*dx #4
        + a_para(chi, mu)*dx
        + a_para(l, tau)*dx #5
        )    
    
    solve(F==0,z,bcs)    
    (p_h, lam_h, q_h, l_h, mu_h) = z.split()
    u_h = Function(V).interpolate(p_h + q_h) 
    u_h.rename("u_h")
    #Save Data
    if save:
        to_save = [u_h]
        SAVE(to_save, f, b, mesh, u_h, dOmegaDVal, Order, "AP", BC_As_Exact)  
    return u_h    

def Solve_PF_Method(f, b, mesh, dOmegaD, dOmegaIN, dOmegaDVal=0,
                    A_para = 1, A_perp = None, Order = 1, eps = 1e-15,
                    save = False, BC_As_Exact = True):
    U = FunctionSpace(mesh, "CG", Order)
    Q = FunctionSpace(mesh, "CG", Order)   
    Z = U*Q
    
    print(f"PF, Z.dim() = {Z.dim()}, Order = {Order}, eps = {eps}.")
    I = Identity(mesh.geometric_dimension())
    if A_perp == None:
        A_perp = I
    
    #Set boundary conditions
    bc0 = DirichletBC(Z.sub(0), dOmegaDVal, dOmegaD)
    bcs = [bc0]
    if (len(dOmegaIN)>0):
        bc2 = DirichletBC(Z.sub(1), 0, dOmegaIN)
        bcs.append(bc2)
    
    #grad_para = lambda alpha: dot(outer(b,b), grad(alpha))
    grad_perp = lambda alpha: dot(I-outer(b,b), grad(alpha))
    a_perp = lambda alpha, beta: inner(A_perp*grad_perp(alpha), grad_perp(beta))
    #a_para = lambda alpha, beta: inner(A_para*grad_para(alpha), grad_para(beta))    
    
    z = Function(Z)
    u, q = split(z)
    v, w = split(TestFunction(Z))

    F = (
        + inner(grad_perp(u), grad_perp(v))*dx
        + inner(q*b, grad(v))*dx
        - inner(f, v)*dx
        - eps * inner(q, w)*dx
        + inner(grad(u), w*b)*dx
        )

    solve(F == 0, z, bcs)  
    (u_h, q_h) = z.split()
    u_h.rename("u_h")
    q_h.rename("q_h")
    if save:
        to_save = [u_h, q_h]
        SAVE(to_save, f, b, mesh, u_h, dOmegaDVal, Order, "PF", BC_As_Exact)
    return u_h

def Solve_PF_STAB_Method(f, b, mesh, dOmegaD, dOmegaIN, dOmegaDVal=0,
                         A_para = 1, A_perp = None, Order = 1, eps = 1e-15,
                         sigma = 0.1, save = False, BC_As_Exact = True):
    U = FunctionSpace(mesh, "CG", Order)
    Q = FunctionSpace(mesh, "CG", Order)   
    Z = U*Q
    
    print(f"PF_STAB, Z.dim() = {Z.dim()}, Order = {Order}, eps = {eps}.")
    I = Identity(mesh.geometric_dimension())
    if A_perp == None:
        A_perp = I
    
    #Set boundary conditions
    bc = DirichletBC(Z.sub(0), dOmegaDVal, dOmegaD)
    
    #grad_para = lambda alpha: dot(outer(b,b), grad(alpha))
    grad_perp = lambda alpha: dot(I-outer(b,b), grad(alpha))
    a_perp = lambda alpha, beta: inner(A_perp*grad_perp(alpha), grad_perp(beta))
    #a_para = lambda alpha, beta: inner(A_para*grad_para(alpha), grad_para(beta))    
    
    z = Function(Z)
    u, q = split(z)
    v, w = split(TestFunction(Z))

    F = (
        + inner(grad_perp(u), grad_perp(v))*dx
        + inner(q*b, grad(v))*dx
        - inner(f, v)*dx
        - eps * inner(q, w)*dx
        - sigma*inner(q, w)*dx
        + inner(grad(u), w*b)*dx
        )

    solve(F == 0, z, bc)  
    (u_h, q_h) = z.split()
    u_h.rename("u_h")
    q_h.rename("q_h")
    if save:
        to_save = [u_h, q_h]
        SAVE(to_save, f, b, mesh, u_h, dOmegaDVal, Order, "PF_STAB", BC_As_Exact) 
    return u_h

def Solve_DN_Method(f, b, mesh, dOmegaD, dOmegaIN, dOmegaDVal=0, A_para = 1, A_perp = None,
                          Order = 1, eps = 1e-15, eps_0 = 1e-7, Repeats = 8):
    V = FunctionSpace(mesh, "CG", Order)
    print(f"DN, Z.dim() = {Z.dim()}, Order = {Order}, eps = {eps}.")
    
    #v, w = TestFunctions(Z)
    v = TestFunction(V)
    w = TestFunction(V)
    I = as_matrix([[1,0],[0,1]])
    
    bc0 = DirichletBC(V, dOmegaDVal, dOmegaD)    
    bcs = [bc0]

    grad_para = lambda u: outer(b,b)*grad(u)
    grad_perp = lambda u: (I-outer(b,b))*grad(u)
    a_para = lambda u, v: inner(grad_para(u),grad(v))*dx 
    a_perp = lambda u, v: inner(grad_perp(u),grad(v))*dx
    a_eps_0 = lambda u, v : a_para(u, v) + eps_0 * a_perp(u, v)

    
    
    u = Function(V).interpolate(x*0)
    q = Function(V).interpolate(x*0)
    
    counter = 0
    L2_errors = []
    H1_errors = []
    while counter < Repeats:
        counter += 1
        u_ = TrialFunction(V) #under score denotes n+1
        solve(a_eps_0(u_, v) == inner(eps_0*f,v)*dx+(eps-eps_0)*a_para(q,v), u, bcs)
        
        q_ = TrialFunction(V)
        solve(a_eps_0(q_, w) == inner(f,w)*dx+eps*a_perp(q,w)-a_perp(u,w) , q, bcs) 
        
        L2_error = norm(u - u_e, "L2")
        H1_error = norm(u - u_e, "H1")
        L2_errors.append(L2_error)
        H1_errors.append(H1_error)
        print("Counter: "+str(counter) + ", ||u - u_exact||_L2: %.1e" % L2_error)
        print("Counter: "+str(counter) + ", ||u - u_exact||_H1: %.1e" % H1_error)    
    return u, L2_errors, H1_errors
