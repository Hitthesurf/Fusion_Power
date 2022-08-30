grad_para = lambda alpha: dot(outer(b,b), grad(alpha))
grad_perp = lambda alpha: dot(I-outer(b,b), grad(alpha))
a_perp = lambda alpha, beta: inner(A_perp*grad_perp(alpha),
                                   grad_perp(beta))
a_para = lambda alpha, beta: inner(A_para*grad_para(alpha), 
                                   grad_para(beta)) 