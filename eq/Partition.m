function [Pi,PiIn,PiOut,ndx] = Partition(ng,p,selNodes,drawThings)
% Partition the sub-graph containing the nodes (selNodes) in a network game (ng) 
% into (numGroups) components. For certain regular structures, the graph is
% divided into approximately equal-sized slices. Otherwise, a spectral
% partitioning algorithm is used.
%
% Outputs:
% Pi (n,numGroups): Pi(i,k) = 1 if node i is in partition k, 0 otherwise
% PiIn (n,numGroups): PiIn(i,k) = 1 if node i is in partition k and has a
%   neighbor in another partition
% PiOut (n,numGroups): PiOut(i,k) = 1 if node i is not in partition k but
%   has a neighbor in partition k
% ndx (1,n): ndx(i) = k means that node i is in partition k

Adj = ng.Adj;
n = size(Adj,1);
nSel = length(selNodes);
numLocalEq = ng.LEQSizes;
markerEdgeWidth = ng.nodeSize/5;

% Partition network
if (strcmp(ng.networkType,'torus') && (ng.nrows == 1 || ng.ncols == 1)) || p == 1
  %  Ring network
  PiGroup = kron(eye(p),ones(ceil(nSel/p),1));
  PiGroup = PiGroup(1:nSel,:);
  [I,J] = find(PiGroup);
  Pi = sparse(selNodes,J,ones(size(I)),n,p);
  [~,ndx] = find(Pi);
elseif strcmp(ng.networkType,'torus') && p == 2
  % Cut in two on longest dimension
  Vx = ng.Vxy(selNodes,1); Vy = ng.Vxy(selNodes,2);
  spanX = max(Vx) - min(Vx);
  spanY = max(Vy) - min(Vy);
  ndx = zeros(n,1);
  if spanX > spanY
    cutX = min(Vx) + (max(Vx) - min(Vx)) / 2;
    ndx(selNodes(Vx<=cutX)) = 1;
    ndx(selNodes(Vx>cutX)) = 2;
  else
    cutY = min(Vy) + (max(Vy) - min(Vy)) / 2;
    ndx(selNodes(Vy<=cutY)) = 1;
    ndx(selNodes(Vy>cutY)) = 2;
  end
  Pi = sparse(selNodes,ndx(selNodes),ones(size(selNodes)),n,2);
elseif strcmp(ng.networkType,'lattice3D') && p == 2
  % Cut in two on longest dimension
  Vx = ng.Vxy(selNodes,1); Vy = ng.Vxy(selNodes,2); Vz = ng.Vxy(selNodes,3);
  spanX = max(Vx) - min(Vx);
  spanY = max(Vy) - min(Vy);
  spanZ = max(Vz) - min(Vz);
  ndx = zeros(n,1);
  [~,longDim] = max([spanX spanY spanZ]);
  switch longDim
    case 1
      cutX = min(Vx) + (max(Vx) - min(Vx)) / 2;
      ndx(selNodes(Vx<=cutX)) = 1;
      ndx(selNodes(Vx>cutX)) = 2;
    case 2
      cutY = min(Vy) + (max(Vy) - min(Vy)) / 2;
      ndx(selNodes(Vy<=cutY)) = 1;
      ndx(selNodes(Vy>cutY)) = 2;
    case 3
      cutZ = min(Vz) + (max(Vz) - min(Vz)) / 2;
      ndx(selNodes(Vz<=cutZ)) = 1;
      ndx(selNodes(Vz>cutZ)) = 2;
  end
  Pi = sparse(selNodes,ndx(selNodes),ones(size(selNodes)),n,2);
  
elseif strcmp(ng.networkType,'tree')
  
  % Fast partitioning algorithm (only for very specific trees!)
  if all(mod(selNodes,2)==1)
    newInds = selNodes + 1;
  else
    newInds = selNodes;
  end
  while all(mod(newInds,2)==0)
    newInds = newInds / 2;
    if all(mod(newInds,2)==1)
      newInds = newInds + 1;
    end
  end
  evenInds = mod(newInds,2)==0;
  oddInds = mod(newInds,2) == 1;
  ndx = zeros(n,1);
  ndx(selNodes(evenInds)) = 1;
  ndx(selNodes(oddInds)) = 2;
  Pi = sparse(selNodes,ndx(selNodes),ones(size(selNodes)),n,2);
  
else % spectral partition
 
  % Symmetrize in case of directed graph
  AdjUndir = (Adj + Adj')/2;
  
  [I,J] = find(AdjUndir);
  
  % Weight graph partition matrix by the numbers of local equilibria
  AdjEq = sparse(I,J,numLocalEq(I)+numLocalEq(J),n,n);
  AdjEq = AdjEq./repmat(sum(AdjEq),n,1);
  subAdj = sparse(AdjEq(selNodes,selNodes));
  
  % Make doubly stochastic
  subAdj=subAdj/((1+eps)*(max(sum(subAdj))));
  subAdj=subAdj+sparse(1:nSel,1:nSel,1-sum(subAdj));
  
  % Partition graph
  [subndx,subPi] = grPartition(subAdj,p,10);
  ndx = zeros(n,1);
  ndx(selNodes) = subndx;
  [I,J] = find(subPi);
  Pi = sparse(selNodes(I),J,subPi(subPi>0),n,p);
end

% Plot partitions
if drawThings
  dim = size(ng.Vxy,2);
  k = max(ndx);
  % colors=hsv(k);
  colors=rand(k,3);
  for i = 1:k
    Vi = find(ndx==i);
    if dim == 2
      plot(ng.Vxy(Vi,1),ng.Vxy(Vi,2),'o','MarkerSize',ng.nodeSize,...
        'LineWidth',markerEdgeWidth,'Color',colors(i,:));
    elseif dim == 3
      plot3(ng.Vxy(Vi,1),ng.Vxy(Vi,2),ng.Vxy(Vi,3),'o','MarkerSize',...
        ng.nodeSize,'LineWidth',markerEdgeWidth,'Color',colors(i,:));
    end
  end
  drawnow;
end

% Allocate group data containers
PiOut = false(n,p); % outer boundary nodes
PiIn = false(n,p); % inner boundary nodes

% For EQ-finding purposes, we use symmetric Adj matrix
AdjSym = ng.Adj | ng.Adj';

for k = 1:p
  
  Vk = find(Pi(:,k));
  VkExp = find(Pi(:,k) | any(AdjSym(:,Vk),2)); 

  VkOut = setdiff(VkExp,Vk);
  PiOut(VkOut,k) = 1;
  
  W = repmat(Pi(:,k),1,length(Vk));
  neighborsK = AdjSym(:,Vk);
  VkInt = Vk(all((neighborsK|W)==W));
  VkIn = setdiff(Vk,VkInt);
  PiIn(VkIn,k) = 1;
  
end

