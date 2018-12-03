function [A,Vxy] = MakeRandGeo(n,p)
% Places vertices randomly on the unit square and connect vertices within a
% distance p of each other

Vxy = rand(n,2);
Dists = dist(Vxy');


A = Dists <= p;
A = A.*(ones(n) - eye(n)); % no self-edges



