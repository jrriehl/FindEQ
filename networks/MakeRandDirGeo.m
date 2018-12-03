function [A,Vxy] = MakeRandDirGeo(n,p)
% Places vertices randomly on the unit square and connect vertices within a
% distance p of each other

Vxy = rand(n,2);
Dists = dist(Vxy');

A = Dists <= p;
A = A - eye(n);

% Directed and no self-edges
[I,J] = find(A);
undir = I < J;
I = I(undir);
J = J(undir);
for x = 1:length(I)
  i = I(x);
  j = J(x);
  if rand > 0.5
    A(i,j) = 0;
  else
    A(j,i) = 0;
  end
end
