#Examples
from firedrake import *
import numpy as np

def Example_1(al = 0., m = 1., Size = 40, eps = 1e-15):
    # Example with no clossed feild lines
    #Create Mesh
    mesh = UnitSquareMesh(Size,Size, quadrilateral= True)
    x, y = SpatialCoordinate(mesh)
    
    #Provide Boundary Conditions
    dOmegaD = [3,4]
    dOmegaIN = [1]
    
    #Calc Vector Field
    b_1 = al*(2*y-1)*cos(m*pi*x)+pi
    b_2 = pi*al*m*(y**2-y)*sin(m*pi*x)
    B = as_vector([b_1,b_2])
    b = B/sqrt(dot(B,B))
    
    #Provide Exact Solution
    u_e0 = sin(pi*y+al*(y**2-y)*cos(m*pi*x))
    u_e1 = cos(2*pi*x)*sin(pi*y)
    u_e = u_e0 + eps*u_e1
    
    #Calc f
    I = Identity(mesh.geometric_dimension())
    A_para = 1
    A_perp = I
    A_eps = (1/eps)*A_para*outer(b,b)
    A_eps += (I-outer(b,b)) * A_perp * (I-outer(b, b))
    f = -div(dot(A_eps, grad(u_e)))
    return f, b, mesh, dOmegaD, dOmegaIN, u_e


def Example_2(eps = 1e-15, Restrict_CL =True):
    # Annulus (Clossed Field Lines)
    # 1<r<2
    mesh = Mesh("Mesh/Annulus_Side_Mesh.msh")
    x, y = SpatialCoordinate(mesh)
    
    #Provide Boundary Conditions
    dOmegaD = [11, 12]
    dOmegaIN = []
    if Restrict_CL:
        dOmegaIN.append(1) #MagSplit Num
    
    #Calc Vector Field
    b_1 = -y
    b_2 = x
    B = as_vector([b_1,b_2])
    b = B/sqrt(dot(B,B))
    
    #Provide Exact Solution
    u_e0 = (x**2+y**2-1)*(4-x**2-y**2)
    u_e1 = (x**2+y**2-1)*(4-x**2-y**2)
    u_e = u_e0 + eps*u_e1
    
    
    #Calc f
    I = Identity(mesh.geometric_dimension())
    A_para = 1
    A_perp = I
    A_eps = (1/eps)*A_para*outer(b,b)
    A_eps += (I-outer(b,b)) * A_perp * (I-outer(b, b))
    f = -div(dot(A_eps, grad(u_e)))
    return f, b, mesh, dOmegaD, dOmegaIN, u_e

def Example_3(a = 0.05, Size = 40,
              eps = 1e-15, Restrict_CL = False):
    #Magnetic Island (Clossed Feild Lines) from DN paper
    #Create Mesh
    mesh = UnitSquareMesh(Size,Size, quadrilateral= True)
    if (Restrict_CL):
        #Load mesh
        mesh = Mesh("Mesh/Square_Mesh.msh")
        
        
    x, y = SpatialCoordinate(mesh)
    
    #Provide Boundary Conditions
    dOmegaD = [3,4]
    dOmegaIN = [1]
    if Restrict_CL:
        dOmegaIN.append(9)
    
    #Calc Vector Field
    b_1 = -cos(pi*y)
    b_2 = 4*a*sin(4*pi*x)
    B = as_vector([b_1,b_2])
    b = B/sqrt(dot(B,B))
    
    #Provide Exact Solution
    u_e0 = sin(10*(sin(pi*y)-a*cos(4*pi*x)))
    u_e1 = cos(2*pi*x)*sin(10*pi*y)
    u_e = u_e0 + eps*u_e1
    
    #Calc f
    I = Identity(mesh.geometric_dimension())
    A_para = 1
    A_perp = I
    A_eps = (1/eps)*A_para*outer(b,b)
    A_eps += (I-outer(b,b)) * A_perp * (I-outer(b, b))
    f = -div(dot(A_eps, grad(u_e)))
    return f, b, mesh, dOmegaD, dOmegaIN, u_e

def Example_4(a = 0.05, Size = 40,
              eps = 1e-15, Restrict_CL = False):
    #Magnetic Island (Clossed Feild Lines) adaption of example 3
    # Enforces vector field to have b dot n = 0, on dOmegaD
    #Create Mesh
    mesh = UnitSquareMesh(Size,Size, quadrilateral= True)
    if (Restrict_CL):
        #Load mesh
        mesh = Mesh("Mesh/Square_Mesh.msh")
        
        
    x, y = SpatialCoordinate(mesh)
    
    #Provide Boundary Conditions
    dOmegaD = [3,4]
    dOmegaIN = [1]
    if Restrict_CL:
        dOmegaIN.append(9)
    
    #Calc Vector Field
    b_1 = -cos(pi*y)
    b_2 = 4*a*sin(4*pi*x)*sin(pi*y)
    B = as_vector([b_1,b_2])
    b = B/sqrt(dot(B,B))
    
    #Provide Exact Solution
    u_e0 = sin(10*sin(pi*y)*e**(-a*cos(4*pi*x)))
    u_e1 = u_e0
    u_e = u_e0 + eps*u_e1
    
    #Calc f
    I = Identity(mesh.geometric_dimension())
    A_para = 1
    A_perp = I
    A_eps = (1/eps)*A_para*outer(b,b)
    A_eps += (I-outer(b,b)) * A_perp * (I-outer(b, b))
    f = -div(dot(A_eps, grad(u_e)))
    return f, b, mesh, dOmegaD, dOmegaIN, u_e 

def Example_5(eps = 1e-15, Restrict_CL = True):
    #Torus (Clossed Field Lines)
    r_I = 1
    mesh = Mesh("Mesh/Torus_Side_Mesh.msh")
    x, y, z = SpatialCoordinate(mesh)
    
    
    #Privide Boundary Conditions
    dOmegaD = [1]
    dOmegaIN = []
    if Restrict_CL:
        dOmegaIN.append(9)
    
    #Calc Vector Field
    b_1 = -y
    b_2 = x
    b_3 = Constant(0)
    B = as_vector([b_1,b_2,b_3])
    b = B/sqrt(dot(B,B))
    
    #Provide Exact Solution
    #r = sqrt(r_I**2+x**2+y**2+z**2-2*r_I*sqrt(x**2+y**2))
    r_s = r_I**2+x**2+y**2+z**2-2*r_I*sqrt(x**2+y**2)
    u_e01 = -r_s + 0.25
    u_e = u_e01 + eps*u_e01
    
    #Calc f
    I = Identity(mesh.geometric_dimension())
    A_para = 1
    A_perp = I
    A_eps = (1/eps)*A_para*outer(b,b)
    A_eps += (I-outer(b,b)) * A_perp * (I-outer(b, b))
    f = -div(dot(A_eps, grad(u_e)))
    return f, b, mesh, dOmegaD, dOmegaIN, u_e
