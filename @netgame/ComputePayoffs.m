function ComputePayoffs(netgame,agents)
% Compute agent payoffs from current state

if nargin < 2
  Edges = netgame.Edges;
else
  edgeInds = netgame.edgeAdjMap(agents,:);
  edgeInds = edgeInds(find(edgeInds));
  Edges = netgame.Edges(edgeInds,:);
end

n = netgame.n;
x = netgame.x;
y = zeros(n,1);
for e = 1:size(Edges,1)
  i = Edges(e,1); j = Edges(e,2);
  y(i) = y(i) + netgame.edgePayoffs(x(i),x(j),e);
end

if netgame.avgPayoffs % Compute average payoffs
  degree = netgame.deg;
  y(degree>0) = y(degree>0)./degree(degree>0);
end

netgame.y = y;

end

