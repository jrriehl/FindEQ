function [eqOpt,compOpt] = ExpandEq(GroupData,targetPhi)
% Recursively expand a tree representation of regional and local 
% equilibrium states into the full set of global equilibrium states.
% Caution: make sure you have enough memory to expand the full set.

tic;
GD = GroupData;
if nargin < 2
  targetPhi = max(GD.uniqPhiPREQ);
end

isTargetPhi = GD.uniqPhiPREQ == targetPhi;
groupInd = find(GD.phiPREQCounts(:,isTargetPhi));
phiInterGroup = GD.phiInterGroup(groupInd);

eqOpt = RecExpandOpt(GD,groupInd,phiInterGroup,targetPhi);
compOpt = toc;

function expEq = RecExpandOpt(GD,groupInd,phiInterGroup,optPhi)
% Return set of all lower level eq that have this score

n = size(GD.Pi,1);

stack = @(x,y) [kron(x,ones(size(y,1),1)) kron(ones(size(x,1),1),y)];
numGroups = size(GD.Pi,2);

if numGroups == 1 || isempty(GD.subGD)
  expEq = GD.PREQ(groupInd,:);
  return;
end

if numGroups > 2
%  error('Not implemented for more than 2 groups yet.');
end

% Otherwise we need to find all possible combinations of the partitioned
% subgraphs whose energy is equal to that of the higher level group
subBREQInd = cell(size(groupInd,1),numGroups);
subBREQ2PREQ = cell(size(groupInd,1),numGroups);
optSubPhi = cell(size(groupInd,1),numGroups);
for i = 1:length(groupInd)
  interPhi = phiInterGroup(i);
  subPhiPREQ = cell(numGroups,1);
  for k = 1:numGroups
    subBREQInd{i,k} = find(GD.PREQ2subBREQ{k}(:,groupInd(i,end)));
    subBREQ2PREQ{i,k} = find(GD.subBREQ2PREQ{k}(:,subBREQInd{i,k}));
    subPhiPREQ{k} = GD.subGD{k}.uniqPhiPREQ;
  end
  [X,Y] = meshgrid(subPhiPREQ{1},subPhiPREQ{2});
  Z = X + Y + interPhi;
  optInd = find(abs(Z-optPhi) < 1e-10);
  if isempty(optInd)
    error('optInd is empty');
  end
  optSubPhi{i,1} = X(optInd);
  optSubPhi{i,2} = Y(optInd);
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
        phiInterGroup = nextGD.phiInterGroup(optGroupInd);        
        nextGD.optInd = stack(groupInd(i,:),optGroupInd);
        eqk = RecExpandOpt(nextGD,optGroupInd,phiInterGroup,optSubPhi{i,k}(j));
      end
      
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
[I,~] = find(GD.Pi);
%WW = setdiff(I,find(sum(expEq,1)));