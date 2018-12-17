function GroupData = JigsawGroups(ng,p,maxGroupSize,selNodes,computeEnergies,drawThings)
% Recursively piece together compatible overlapping local level equilibria
% (smaller groups of nodes) to generate broader equilibria (larger groups
% of nodes) by partitioning 
%
% Inputs
% ng: network game object
% p: number of groups in each partition (usually 2)
% maxGroupSize: maximimum number of nodes in a base-level partition (groups
%   of this size or smaller will not be partitioned further)
% selNodes: defines subgraph corrsesponding to the active group in the recursion
% computeEnergies: if true, keep track of the energy values for the
%   regional equilibrium states and the relationship between groups and
%   partition levels
%
% Output
% GroupData: collection of all relevant information about the results for
%   active group, and lower level groups (allGD is a cell-array of GroupData
%   results from lower levels of partitioning or base level).

% If not specified, include all nodes (means this is probably the top level)
if nargin < 3
  selNodes = 1:size(ng.eqLocal,1);
end

% If the number of nodes is greater than the max group size, then partition the graph
if length(selNodes) > maxGroupSize
  
  % Pi (n,numGroups): Pi(i,k) = 1 if node i is in partition k, 0 otherwise
  % PiIn (n,numGroups): PiIn(i,k) = 1 if node i is in partition k and has a
  %   neighbor in another partition
  % PiOut (n,numGroups): PiOut(i,k) = 1 if node i is not in partition k but
  %   has a neighbor in partition k
  [Pi,PiIn,PiOut] = Partition(ng,p,selNodes,drawThings);
  
  % Allocate subgroup data containers
  subPREQ = cell(p,1);            % PREQ (partition regional eq. states) for subgroups (unique eq. states defined on nodes that lie on a boundary of the lower-level partition)
  subBREQ = cell(p,1);            % BREQ (boundary regional eq. states) for subgroups (unique eq. states defined on nodes that lie on a boundary of the current level partition)
  indSubPREQ2subBREQ = cell(p,1);    % index mapping from subPREQ to subBREQ
  subGD =  cell(p,1);             % GroupData results from partitioned subgroups
  subBREQCountsTotal = cell(p,1); % total number of ways each subBREQ state can be instantiated by lower level configurations
  subBREQ2subPREQ = cell(p,1);       % mapping from subBREQ to (current level) PREQ states
  subUniqPhiBREQ = cell(p,1);     % unique energy values of subBREQ elements
  subPhiBREQCounts = cell(p,1);   % how many subBREQ elements have each unique energy
  subPhiBREQOffset = zeros(p,1);  % offset to make all energy values positive for each subgroup
  
  for k = 1:p
    
    Vk = find(Pi(:,k));
    VkBound = find(PiIn(:,k) | PiOut(:,k));
    Vk = SortNodes(ng,Vk);
    
    % Recursively decompose this group into smaller pieces and solve
    nGroup = min(p,length(Vk));
    if ng.showStatus, disp(['Jigsaw Group: ' num2str(length(Vk))]); end
    subGD{k} = JigsawGroups(ng,nGroup,maxGroupSize,Vk,computeEnergies,drawThings);
    VkPart = find(any(subGD{k}.PiIn + subGD{k}.PiOut,2));
    
    % Reduce dimension by storing only boundary states, while keeping track
    % of how many ways they can be put together
    if ng.showStatus, disp('Mapping to boundary states...'); end
    subPREQ{k} = subGD{k}.PREQ(:,ismember(VkPart,VkBound));
    [subBREQ{k},X,indSubPREQ2subBREQ{k}] = unique(subPREQ{k},'rows');
    
    % While reducing to boundary states, keep track of the number of group
    % EQ configurations included by each configuration of boundary EQ
    I = (1:length(indSubPREQ2subBREQ{k}))';
    subBREQ2subPREQ{k} = sparse(I,indSubPREQ2subBREQ{k},subGD{k}.PREQCountsTotal,length(indSubPREQ2subBREQ{k}),length(X));
    subBREQCountsTotal{k} = sum(subBREQ2subPREQ{k})';
    subBREQ2subPREQ{k} = subBREQ2subPREQ{k} > 0;
    
    if computeEnergies && ~isempty(subGD{k}.uniqPhiPREQ)
      % Shift range of energies to be positive to avoid important zero entries
      subPhiBREQOffset(k) = -min(subGD{k}.uniqPhiPREQ) + 1;
      
      % subPhiBREQCounts{k}(i,j) = number of ways the unique energy value
      % subGD{k}.uniqPhiPREQ(i) is achieved within the jth subBREQ{k} entry
      pPCk = subGD{k}.phiPREQCounts';
      subPhiBREQCounts{k} = zeros(length(subGD{k}.uniqPhiPREQ),length(X));
      for i = 1:length(indSubPREQ2subBREQ{k})
        subPhiBREQCounts{k}(:,indSubPREQ2subBREQ{k}(i)) = subPhiBREQCounts{k}(:,indSubPREQ2subBREQ{k}(i)) + pPCk(:,i);
      end
      subPhiBREQCounts{k} = sparse(subPhiBREQCounts{k})';
         
      % subUniqPhiBREQ{k} = vector of all unique feasible energy values for partition group k
      subUniqPhiBREQ{k} = subGD{k}.uniqPhiPREQ;
    end
    
  end
  if ng.showStatus, disp('Group equilibria computed.'); end
  
  % Partition matrix including inter-group boundary nodes at current level
  PiPart = PiIn | PiOut;
  
  % Since we typically use only two groups, this order does not matter much.
  % When partitioning into more than two groups, it may become more important.
  groupOrder = 1:size(subBREQ,1);
  
  % Compute intersections of boundary regional equilibrium sets subEqBound
  % PREQ (nPREQ,nVpart): partion reg. eq. states
  % Vpart: nodes that lie on a boundary of the current level partition
  if ng.showStatus, disp(['Jigsaw Meta: ' num2str(length(groupOrder))]); disp(subBREQ); end
  PREQ = GroupMerge(subBREQ,groupOrder,PiPart);
  Vpart = find(any(PiIn+PiOut,2));
  if ng.showStatus, disp('Group connections computed.'); end
  
  % Extra trimming step for 3D lattice graphs (TODO: generalize)
  Vout = setdiff(Vpart,selNodes);
  if strcmp(ng.networkType,'lattice3D') && ~isempty(Vout) && ~isempty(PREQ)
    PREQ = TrimRegionalEqLattice3D(ng,PREQ,Vpart,Vout);
  end
  
  % phiInterGroup (nPREQ,1): total energy of edges that connect different 
  % subgroups for each PREQ state
  if ng.showStatus, disp(['Computing inter-group energies: ' num2str(size(PREQ,1))]); end
  if computeEnergies
    if isempty(PREQ)
      phiInterGroup = [];
    else
      if strcmp(ng.updateRules{1},'BestResponseControl')
        phiInterGroup = sparse(size(PREQ,1),1);
      else
        phiInterGroup = ComputeGroupEnergies(ng,ng.Adj,PREQ,Pi,Vpart);
      end
    end
  end
  
  if ng.showStatus
    disp(['Calculating number of equilibria: ' num2str(size(PREQ,1)*p)]); 
  end
  
  % Divide PREQ and BREQ into full matrices of group components (for processing speed)
  PREQDiv = cell(p,1);
  subBREQDiv = cell(p,1);
  for k = 1:p
    Vk = find(PiPart(:,k));
    PREQDiv{k} = PREQ(:,ismember(Vpart,Vk));
    subBREQk = uint8(full(subBREQ{k}(:,ismember(Vk,Vpart))));
    subBREQDiv{k} = reshape(subBREQk.',1,length(Vk),[]);
  end
  
  % Compute the total number of equilibria without actually computing them all
  % indPREQ2subBREQ: index mapping from PREQ to subBREQ states in each subgroup
  % PREQCountsTotal: total number of ways that each PREQ state can be instantiated by lower-level states
  % PREQ2subBREQ {p,1}(nPREQ,nSubBREQ): PREQ2subBREQ(i,j) = 1 if nPREQ state i is consistent with subBREQ state j  
  indPREQ2subBREQ = zeros(size(PREQ,1),p);
  PREQCountsTotal = zeros(size(PREQ,1),p);
  PREQ2subBREQ = cell(p,1);
  
  for k = 1:p
    PREQ2subBREQ{k} = AllBsxSubDivParallel(@eq,subBREQDiv{k},PREQDiv{k})';
    [indPREQ2subBREQ(:,k),~] = find(PREQ2subBREQ{k});
    if isempty(indPREQ2subBREQ(:,k))
      PREQCountsTotal(:,k) = 0;
      uniqPhiPREQ = zeros(0,0);
      phiPREQCounts = zeros(0,0);
    else
      PREQCountsTotal(:,k) = subBREQCountsTotal{k}(indPREQ2subBREQ(:,k));
    end
  end
  numEq = sum(prod(PREQCountsTotal,2));
  if computeEnergies && ~isempty(subUniqPhiBREQ)
    [uniqPhiPREQ,phiPREQCounts] = MapGroupEnergies(subUniqPhiBREQ,subPhiBREQCounts,indPREQ2subBREQ,phiInterGroup,subPhiBREQOffset);
  end
  
  if ng.showStatus, disp('Done'); end
  
  % Compute number of edges but by partition
  VinQ = find(any(PiIn,2)); 
  VoutQ = find(any(PiOut,2));
  numCutEdges = sum(ismember(ng.Edges(:,1),VinQ) & ismember(ng.Edges(:,2),VoutQ)); 
  
  % Return partition data
  GroupData.Pi = Pi;
  GroupData.PiOut = PiOut;
  GroupData.PiIn = PiIn;
  GroupData.numCutEdges = numCutEdges;
  GroupData.PREQ = PREQ;
  GroupData.subBREQCountsTotal = subBREQCountsTotal;
  GroupData.subGD = subGD;
  GroupData.PREQCountsTotal = prod(PREQCountsTotal,2);
  GroupData.numEq = numEq;
  if computeEnergies
    GroupData.subUniqPhiBREQ = subUniqPhiBREQ;
    GroupData.subPhiBREQCounts = subPhiBREQCounts;
    GroupData.uniqPhiPREQ = uniqPhiPREQ;
    GroupData.phiPREQCounts = phiPREQCounts;
    GroupData.phiInterGroup = phiInterGroup;
    GroupData.subBREQ2PREQ = subBREQ2subPREQ;
    GroupData.PREQ2subBREQ = PREQ2subBREQ;
  end
  
else % base level - just compute all compatible equilibria within this group without partitioning further
  
  nGroup = length(selNodes);
  if ng.showStatus, display(['Jigsaw Local: ' num2str(nGroup)]); end
  
  selNodes = SortNodes(ng,selNodes);
  
  % Construct partition matrix as if each node is it's own partition group
  Pi = false(ng.n,nGroup);
  PiPart = Pi;
  Pi(selNodes,1:nGroup) = eye(nGroup,nGroup);
  AdjSym = ng.Adj | ng.Adj';
  for ii = 1:length(selNodes)
    i = selNodes(ii);
    PiPart(:,ii) = Pi(:,ii) | AdjSym(:,i)>0;
  end
  
  if ng.showStatus, display(ng.LEQ(selNodes)); end
  [PREQ,Vpart] = GroupMerge(ng.LEQ,selNodes,PiPart);
  Vpart = find(Vpart);
  
  % TEMP for proteins
  Vout = setdiff(Vpart,selNodes);
  if ~isempty(Vout) && strcmp(ng.networkType,'lattice3D')
    PREQ = TrimRegionalEqLattice3D(ng,PREQ,Vpart,Vout);
  end
  
  nPREQ = size(PREQ,1);
  GroupData.Pi = Pi;
  GroupData.PiOut = PiPart;
  GroupData.PiIn = PiPart;
  GroupData.numCutEdges = 0;
  GroupData.PREQ = PREQ;
  GroupData.subPhiBREQCounts = repmat({1},length(selNodes),1);
  GroupData.subGD = [];
  GroupData.PREQCountsTotal = ones(nPREQ,1);
  GroupData.numEq = nPREQ;
  if computeEnergies
    phiPREQ = ComputeEnergy(ng,ng.Adj,PREQ,Vpart,selNodes);
    GroupData.phiInterGroup = zeros(nPREQ,1);
    GroupData.uniqPhiPREQ = unique(phiPREQ);
    phiPREQCounts = zeros(nPREQ,size(GroupData.uniqPhiPREQ,1));
    for i  = 1:length(GroupData.uniqPhiPREQ)
      phiPREQCounts(:,i) = phiPREQ == GroupData.uniqPhiPREQ(i);
    end
    GroupData.phiPREQCounts = phiPREQCounts;
    GroupData.subBREQ2PREQ = {sparse(1:nPREQ,1:nPREQ,ones(nPREQ,1))};
    GroupData.PREQ2subBREQ = {sparse(1:nPREQ,1:nPREQ,ones(nPREQ,1))};
  end
  
end

end

function sorted = SortNodes(ng,nodeList)
% The first node is the one having the fewest local EQ. The rest are sorted
% by distance to this node (primary index) and number of LEQ (secondary).

numLEQ = ng.LEQSizes(nodeList);
[~,minS] = min(numLEQ);
dist = bfs(ng.Adjs,nodeList(minS));
dist2root = dist(nodeList);
dist2root(dist2root==-1) = inf;
[~,order] = sortrows([dist2root -numLEQ]);
sorted = nodeList(order);

end


