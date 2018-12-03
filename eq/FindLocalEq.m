function LEQ = FindLocalEq(ng)
% Return a cell array containing the set of all local equilibrium states in
% the respective neighborhoods of each node in the network. The resulting
% sets are indexed by the order they appear columnwise in the adj matrix.

Adj = ng.Adj;
Adjs = sparse(Adj);
n = size(Adj,1);
allS = (1:ng.ns)';
% AdjN = sparse(Adj | Adj' + eye(n));

% First, compute all locally stable strategy states
LEQ = cell(n,1);
neighbors1 = cell(n,1);
neighbors2 = cell(n,1);
numEqs = zeros(n,1);

for i = 1:n
  
  ngi = ng;
  
  if ng.showStatus
    display(['Computing local eq. at node: ' num2str(i)]);
  end
   
  % Find all (in and out) neighbors
  neighbors = find(Adjs(:,i) + Adjs(i,:)');
  neighbors1{i} = [i;neighbors];  
  
  % Neighborhood depends on update rule
  switch ng.updateRules{i}
    case {'Imitation','FastIM'}
      J = find(any(Adjs(neighbors,:) + Adjs(:,neighbors)',1));      
      nextNeighbors = setdiff(unique(J),i)';
      neighbors2{i} = [i;unique([neighbors;nextNeighbors])];
      neighborhood = neighbors2{i};
    otherwise
      neighborhood = neighbors1{i};
  end
  numNeighbors = length(neighborhood);
  
  % Generate matrix of all possible local state configurations
  Li = LogMat(numNeighbors,allS);
  
  % Find all locally invariant neighbor configurations
  isStable = false(size(Li,1),1);
  for j = 1:size(Li,1)
    x = ones(n,1);
    if ismember(ng.updateRules{i},{'Imitation','FastIM'})
      x(neighbors2{i}) = Li(j,:);
      ngi.x = x;
      ComputePayoffs(ng,neighbors2{i});
    else
      x(neighbors1{i}) = Li(j,:);
      ngi.x = x;
    end
    Update(ngi,i);
    xplus = ngi.x;
    if x(i) == xplus(i)
      isStable(j) = true;
    end
  end
  Ls = Li(isStable,:);
  numEqs(i) = sum(isStable);
  [~,order] = sort(neighborhood);
  LEQ{i} = Ls(:,order);  
end

ng.LEQ = LEQ;
ng.LEQSizes = arrayfun(@(x)size(LEQ{x},1),(1:n)');