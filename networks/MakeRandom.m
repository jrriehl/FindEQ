function [A,Vxy] = MakeRandom(n,p)
% Generate Erdos-Renyi random graph
% Vxy = rand(n,2);
theta = 2*pi*(1:n)'/n;
Vxy = [cos(theta) sin(theta)];

A = rand(n,n) <= p;
A = A.*(ones(n) - eye(n)); % no self-edges

[I,J] = find(ones(n));
A(I>J) = 0;
A = A + A';



