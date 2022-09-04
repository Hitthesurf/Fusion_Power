V = FunctionSpace(mesh, "CG", Order)
...
bcs = [DirichletBC(V, dOmegaDVal, dOmegaD)] 
...
u = TrialFunction(V)
psi = TestFunction(V)
u_h = Function(V)
a = ((1/eps)*a_para(u,psi) + a_perp(u, psi))*dx
F = inner(f, psi)*dx
solve(a==F, u_h, bcs)