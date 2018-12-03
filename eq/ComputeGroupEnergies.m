function out = ComputeGroupEnergies(ng,Adj,PREQ,Pi,Vpart)
% Compute cost function only on edges that connect nodes in different
% groups of the partition Pi. Divide into smaller pieces for parallel
% computation.

numEq = size(PREQ,1);
numGroups = size(Pi,2);
AdjOut = Adj;
maxSize = 1e7;

for k = 1:numGroups
  Vk = find(Pi(:,k));
  AdjOut(Vk,Vk) = 0;
end

nonSelNodes = sum(Pi,2)==0;
AdjOut(nonSelNodes,:) = 0;
AdjOut(:,nonSelNodes) = 0;

[neq,nv] = size(PREQ);

if neq*nv < maxSize % too small, don't bother to subdivide
  out = ComputeEnergy(ng,AdjOut,PREQ,Vpart,Vpart);
else
  
  nDiv = ceil(neq*nv / maxSize);
  divInds = round(linspace(1,numEq+1,nDiv+1));
  
  PREQSlices = cell(nDiv,1);
  for k = 1:nDiv
    divIndsk = divInds(k):divInds(k+1)-1;
    PREQSlices{k} = PREQ(divIndsk,:);
  end
 
  phiAll = cell(nDiv,1);
  parfor k = 1:nDiv
    display(['GroupEnergy: ' num2str(k) ' of ' num2str(nDiv)]);
    phiAll{k} = ComputeEnergy(ng,AdjOut,PREQSlices{k},Vpart,Vpart);
  end
  
  out = cat(1,phiAll{:});
  
end