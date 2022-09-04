bcs = [bc0]
...
F = (+inner(grad_perp(u), grad_perp(v))*dx+inner(q*b, grad(v))*dx
     - inner(f, v)*dx #1
     -eps*inner(q, w)*dx-sigma*inner(q, w)*dx
     + inner(grad(u), w*b)*dx) #2