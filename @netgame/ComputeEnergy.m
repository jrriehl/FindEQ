function out = ComputeEnergy(netgame,Adj,X,Xind,selNodes)

n = size(Adj,1);
if nargin < 4
  selNodes = (1:n)';
end

selNodes = sort(selNodes);

out = netgame.energyFunction(Adj,X,Xind,selNodes);