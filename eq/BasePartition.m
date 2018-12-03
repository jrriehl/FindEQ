function [Pi,PiPart] = BasePartition(Adj,selInds)
% Construct partition matrix as if each node is it's own partition group

n = size(Adj,1);
nGroup = length(selInds);

% In case of directed graphs, we need to include both in and out neighbors
Adj = Adj | Adj';

Pi = false(n,nGroup);
PiPart = Pi;
Pi(selInds,1:nGroup) = eye(nGroup,nGroup);
for ii = 1:length(selInds)
  i = selInds(ii);
  PiPart(:,ii) = Pi(:,ii) | (Adj(i,:)>0)';
end