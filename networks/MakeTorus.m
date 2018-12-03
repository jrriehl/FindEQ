function [A,Vxy] = MakeTorus(nrows,ncols)
% Make a lattice that wraps around itself - a torus

n = nrows*ncols;
V = (1:n)';
A = zeros(n,n);

Edown = [V mod(V+ncols-1,n)+1];
Eright = [V ncols*floor((V-1)/ncols) + mod(V,ncols)+1];

E = [Edown;Eright];

for e = 1:size(E,1)
  i = E(e,1);
  j = E(e,2);
  if i==j
    continue; % No self edges
  end
  A(i,j) = 1;
  A(j,i) = 1;
end

qr = (1/nrows:1/nrows:1) - 1/(2*nrows);
qc = (1/ncols:1/ncols:1) - 1/(2*ncols);
Vx = kron(qr,ones(ncols,1));
Vy = kron(qc,ones(nrows,1))';
Vxy = [Vx(:) Vy(:)];

