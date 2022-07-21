#Euro_Fusion_Example_1
from firedrake import *
import numpy as np
import sympy as sp


def Get_Func_Repeat(eps = 1e-15, eps_0 = 0.001, Order = 1, al = 0., m = 1., Repeats = 4, Size = 40):
    
    #mesh = Mesh('Annulus_Mesh005.msh')
    mesh = UnitSquareMesh(Size,Size, quadrilateral= True)
    V = FunctionSpace(mesh, "CG", Order)
    W = FunctionSpace(mesh, "CG", Order)
    Z = V*W

    
    #v, w = TestFunctions(Z)
    v = TestFunction(V)
    w = TestFunction(W)
    I = as_matrix([[1,0],[0,1]])
    x, y = SpatialCoordinate(mesh)
    
    #Now Exact Function
    u_e = sin(pi*y+al*(y**2-y)*cos(m*pi*x)) + eps*cos(2*pi*x)*sin(pi*y)

    B = as_vector([al*(2*y-1)*cos(m*pi*x)+pi,pi*al*m*(y**2-y)*sin(m*pi*x)])
    b = B/sqrt(dot(B,B))

    
    #Calc f
    para_u = div(outer(b,b)*grad(u_e))
    perp_u = div((I-outer(b,b))*grad(u_e))
    f = -para_u/eps -perp_u
    
    bc0 = DirichletBC(V, 0, 3)
    bc1 = DirichletBC(V, 0, 4)
    
    bcs = [bc0, bc1]

    grad_para = lambda u: outer(b,b)*grad(u)
    grad_perp = lambda u: (I-outer(b,b))*grad(u)
    a_para = lambda u, v: inner(grad_para(u),grad(v))*dx 
    a_perp = lambda u, v: inner(grad_perp(u),grad(v))*dx
    a_eps_0 = lambda u, v : a_para(u, v) + eps_0 * a_perp(u, v)

    
    
    u = Function(V).interpolate(x*0)
    q = Function(W).interpolate(x*0)
    
    counter = 0
    while counter < Repeats:
        counter += 1
        u_ = TrialFunction(V) #under score denotes n+1
        solve(a_eps_0(u_, v) == inner(eps_0*f,v)*dx+(eps-eps_0)*a_para(q,v), u, bcs)
        
        q_ = TrialFunction(W)
        solve(a_eps_0(q_, w) == inner(f,w)*dx+eps_0*a_perp(q,w)-a_perp(u,w) , q) 

      
        L2_error = norm(u - u_e, "L2")
        H1_error = norm(u - u_e, "H1")
        print("Counter: "+str(counter) + ", ||u - u_exact||_L2: %.1e" % L2_error)
        #print("||u - u_exact||_H1: %.1e" % H1_error)
    
    return u, q, u_e