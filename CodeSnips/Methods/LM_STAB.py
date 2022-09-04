V = FunctionSpace(mesh, "CG", Order) 
Z = V*V
...
bc0 = DirichletBC(Z.sub(0), 0, dOmegaD)
bc1 = DirichletBC(Z.sub(1), 0, dOmegaD)
bcs = [bc0, bc1]
...
z = Function(Z)
u, q = split(z) #Trial Function
u_, q_ = split(TestFunction(Z))
F = (+ a_perp(u, u_)*dx+a_para(q, u_)*dx-inner(f, u_)*dx
     + a_para(u, q_)*dx-sigma*inner(q,q_)*dx) 
solve(F==0,z,bcs)    
(u_h, q_h) = z.split()