function [E,c] = computeBands(V,j,q);

n=length(j);
Vmat = (V/2)*(sparse(1:n,1:n,1)+.5*sparse(2:n,1:n-1,1,n,n)+.5*sparse(1:n-1,2:n,1,n,n));
Tmat = sparse(diag((q-2*j).^2));
H = Tmat+Vmat;
[c,E]=eig(full(H));
E = diag(E);