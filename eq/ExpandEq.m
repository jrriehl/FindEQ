function [eqOpt,compOpt] = ExpandEq2(GroupData,targetPhi)
% Recursively expand a tree representation of regional and local equilibrium 
% states into the set of global equilibrium states whose cost is targetPhi.
% Caution: make sure you have enough memory to expand the full set.

tic;
GD = GroupData;
if nargin < 2
  targetPhi = max(GD.uniqPhiPREQ);
end

isTargetPhi = GD.uniqPhiPREQ == targetPhi;

% Initialize the algorithm with the indices of top-level PREQ sets that can
% achieve the cost targetPhi
groupInd = find(GD.phiPREQCounts(:,isTargetPhi));

% Extract the inter-partition costs for these indices
phiInterGroup = GD.phiInterGroup(groupInd);

% Initiate the recursion to expand all eq. states whose cost is targetPhi
eqOpt = RecExpandOpt(GD,groupInd,phiInterGroup,targetPhi);
compOpt = toc;

function expEq = RecExpandOpt(GD,groupInd,phiInterGroup,optPhi)
% Return set of all lower level eq that have this score

n = size(GD.Pi,1);
numGroups = size(GD.Pi,2);

% Function that generates rows consisting of all combinations of the rows of x and y
stack = @(x,y) [kron(x,ones(size(y,1),1)) kron(ones(size(x,1),1),y)];

% If this is the base level, simply return the PREQ set which is the full REQ set
if numGroups == 1 || isempty(GD.subGD)
  expEq = GD.PREQ(groupInd,:);
  return;
end

% Otherwise we need to find all possible combinations of the partitioned
% subgraphs whose cost is equal to that of the higher level group
subBREQ2PREQ = cell(size(groupInd,1),numGroups);
optSubPhi = cell(size(groupInd,1),numGroups);
for i = 1:length(groupInd)
  
  interPhi = phiInterGroup(i);
  subPhiPREQ = cell(numGroups,1);
  for k = 1:numGroups
    % Get the corresponding indices of the subBREQ states
    subBREQInd = GD.PREQ2subBREQ{k}(:,groupInd(i,end));
    % Get the corresponding indices of the subPREQ states
    subBREQ2PREQ{i,k} = find(GD.subBREQ2PREQ{k}(:,subBREQInd));
    % Get the unique cost values for the subPREQ states
    subPhiPREQ{k} = GD.subGD{k}.uniqPhiPREQ;
  end
  
  % Construct a matrix whose entries represent the total costs of all
  % combinations of subPREQ states. This will allow us to recover the
  % indices of the subPREQ states that can attain the desired cost.
  expSubPhiPREQ = repmat({ones(size(subPhiPREQ{1},1),1)},numGroups,1);
  expSubPhiPREQ{1} = subPhiPREQ{1};
  numSubPhis = arrayfun(@(x)size(subPhiPREQ{x},1),(1:numGroups)');
  expPhiPREQ = zeros(prod(numSubPhis),1);  
  for k = 1:numGroups
    for l = 2:numGroups
      if k==l
        expSubPhiPREQ{k} = kron(expSubPhiPREQ{k},subPhiPREQ{l});
      else
        expSubPhiPREQ{k} = kron(expSubPhiPREQ{k},ones(size(subPhiPREQ{l},1),1));
      end
    end
    expPhiPREQ = expPhiPREQ + expSubPhiPREQ{k};
  end
  expPhiPREQ = expPhiPREQ + interPhi;

  % Find the indices of the combinations that yield the desired cost
  optInd = abs(expPhiPREQ-optPhi) < 1e-10;
  
  % Get the costs of the subPREQ states corresponding to these combinations
  for k = 1:numGroups
    optSubPhi{i,k} = expSubPhiPREQ{k}(optInd);
  end
end

allOptEq = [];
for i = 1:length(groupInd)
  
  numPhiConfigs = length(optSubPhi{i,1});
  
  for j = 1:numPhiConfigs
    
    groupEq = zeros(0,1);
    Vgroup = false(n,1);
    skip = false;
    
    for k = 1:numGroups
      
      nextGD = GD.subGD{k};
      ViExp = GD.Pi(:,k) | GD.PiOut(:,k);
      
      % Select indices of partitioned EQ states that (i) are
      % compatible with the higher level BEQ state and (ii) match the
      % specified cost value optSubPhi{i,k}(j)
      optUniqPhiInd = ismember(nextGD.uniqPhiPREQ,optSubPhi{i,k}(j));
      candOptSubPhiInd = nextGD.phiPREQCounts(subBREQ2PREQ{i,k},optUniqPhiInd) > 0;
      optGroupInd = subBREQ2PREQ{i,k}(candOptSubPhiInd);
      
      % In some cases, this particular configuration may not be feasible
      if isempty(optGroupInd)
        skip = true;
        break;
      else
        % Otherwise, continue the recursive expansion until we reach the base level
        phiInterGroup = nextGD.phiInterGroup(optGroupInd);        
        nextGD.optInd = stack(groupInd(i,:),optGroupInd);
        eqk = RecExpandOpt(nextGD,optGroupInd,phiInterGroup,optSubPhi{i,k}(j));
      end
      
      % Merge the resulting REQ sets
      if isempty(groupEq)
        groupEq = eqk;
      else
        groupEq = Merge(groupEq,eqk,Vgroup,ViExp);
      end
      Vgroup = Vgroup | ViExp;
    end
    
    if ~skip
      allOptEq = [allOptEq;groupEq];
    end
  end
end

expEq = allOptEq;