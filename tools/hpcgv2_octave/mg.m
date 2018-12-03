## A three level multigrid solver demonstrating main aspects of HPCG

function [z] =  mg(r0,sol,...
                         A0,A1,A2,A3,...
                         L0,L1,L2,L3,...
                         D0,D1,D2,D3,...
                         U0,U1,U2,U3,n)
                         
n1=n/2; n2=n1/2; n3=n2/2;
%% Multigrid
s0 = cg(r0, L0, sol);
s0=D0*s0;
z0=CG(r0, U0, s0);
%%%%
z0=reshape(z0,n,n,n);
sol1(1:n1,1:n1,1:n1)=z0(1:2:n,1:2:n,1:2:n);
sol1=reshape(sol1,n1^3,1);
s1=cg(r1, L1, sol1);
s1=D1*s1;
z1=cg(r1, U1, s1);
%%%%
sol2(1:n2,1:n2,1:n2)=z1(1:2:n1,1:2:n1,1:2:n1);
sol2=reshape(sol2,n2^3,1);
s2=cg(r2, L2, sol2);
s2=D2*s2;
z2=cg(r2, U2, s2);
%%%%
sol3(1:n3,1:n3,1:n3)=z2(1:2:n2,1:2:n2,1:2:n2);
sol3=reshape(sol3,n3^3,1);
s3=cg(r3, L3, sol3);
s3=D3*s3;
z3=cg(r3, U3, s3);
%%%%
%%%%
z3=reshape(z3,n3,n3,n3);
sol2=reshape(sol2,n2,n2,n2);
sol2(1:2:n2,1:2:n2,1:2:n2)=sol2(1:2:n2,1:2:n2,1:2:n2)...
                          +z3(1:n3,1:n3,1:n3);
sol2=reshape(sol2,n2^3,1);
s2=cg(r2, L2, sol2);
s2=D2*s2;
z2=cg(r2, U2, s2);
%%%%
z2=reshape(z2,n2,n2,n2);
sol1=reshape(sol1,n1,n1,n1);
sol1(1:2:n1,1:2:n1,1:2:n1)=sol1(1:2:n1,1:2:n1,1:2:n1)...
                          +z2(1:n2,1:n2,1:n2);
sol1=reshape(sol1,n1^3,1);
s1=cg(r1, L1, sol1);
s1=D1*s1;
z1=cg(r1, U1, s1);
%%%%
z1=reshape(z1,n1,n1,n1);
sol0=reshape(sol0,n,n,n);
sol0(1:2:n,1:2:n,1:2:n)=sol0(1:2:n,1:2:n,1:2:n)...
                       +z1(1:n1,1:n1,1:n1);
sol0=reshape(sol0,n^3,1);
s0=cg(r0, L0, sol0);
s0=D0*s0;
z=cg(r0, U0, s0);
endfunction
