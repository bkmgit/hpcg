% An example program to demonstrate main operations in HPCG with multigrid preconditioner
function [sol] = hpcg(f0,...
                             A0,A1,A2,A3,...
                             L0,L1,L2,L3,...
                             D0,D1,D2,D3,...
                             U0,U1,U2,U3,n)


% initial guess
sol=zeros(size(f));
% initial residual 
r=f0-A0*sol; 
% do first iteration outside of loop to avoid if 
% statement

z=mg(r,sol,...
           A0,A1,A2,A3,...
           L0,L1,L2,L3,...
           D0,D1,D2,D3,...
           U0,U1,U2,U3,n);

%%
p=z;
alpha=dot(r,z);
ap=A*p;
alpha=alpha/dot(p,ap);
sol=sol+alpha*p;
r=r-alpha*ap;
for k=0:20
  %% Multigrid
  z=mg(r,sol,...
           A0,A1,A2,A3,...
           L0,L1,L2,L3,...
           D0,D1,D2,D3,...
           U0,U1,U2,U3,n);
  %%%
  p=z;
  alphaold=alpha;
  alpha=dot(r,z);
  beta=alpha/alphaold;
  p=beta*p+z;
  ap=A*p;
  alpha=alpha/dot(p,ap);
  sol=sol+alpha*p;
  r=r-alpha*ap;
  if (norm(r) < 1e-6)|(alpha < 1e-6), break, end;
endfor

endfunction
