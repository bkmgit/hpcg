function [sol] = CG(r, A, sol)
p=r;
% Perform a single CG iteration
beta=dot(r,r);
ap=A*p;
alpha=beta/dot(p,ap);
sol=sol+alpha*p;

endfunction
