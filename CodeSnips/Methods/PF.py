V = FunctionSpace(mesh, "CG", Order)
Q = FunctionSpace(mesh, "CG", Order)   
Z = V*Q
...
bc0 = DirichletBC(Z.sub(0), 0, dOmegaD)
bc2 = DirichletBC(Z.sub(1), 0, dOmegaIN)
bcs = [bc0, bc2]
...
z = Function(Z)
u, q = split(z)
v, w = split(TestFunction(Z))
F = (+inner(grad_perp(u), grad_perp(v))*dx+inner(q*b, grad(v))*dx
     - inner(f, v)*dx #1
     - eps * inner(q, w)*dx + inner(grad(u), w*b)*dx) #2
solve(F == 0, z, bcs)  
(u_h, q_h) = z.split()