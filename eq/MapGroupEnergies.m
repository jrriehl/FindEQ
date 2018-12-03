function [uniqPhiPREQ,phiPREQCounts] = MapGroupEnergies(subUniqPhiBREQ,...
  subPhiBREQCounts,indPREQ2subBREQ,phiInterGroup,subPhiBREQOffset)
% Generate the matrices that map the energy values between levels in the decomposition

maxSize = 1e6;
numGroups = length(subUniqPhiBREQ);

subPhiPREQ = cell(numGroups,1);
subPhiPREQCounts = cell(numGroups,1);
for k = 1:numGroups
  % subPhiPREQ{k}(i,j) = energy value uniqPhiBound{k}(i) if achievable
  % within lower-level eqBound entry j (zero otherwise)
  uniqPhiBoundShift = subUniqPhiBREQ{k} + subPhiBREQOffset(k);
  Z = (subPhiBREQCounts{k}>0).*repmat(uniqPhiBoundShift',size(subPhiBREQCounts{k},1),1);
  subPhiPREQ{k} = Z(indPREQ2subBREQ(:,k),:);
  
  % subPhiPREQCounts{k}(i,j) = number of ways the unique energy value
  % uniqPhiBound{k}(i) is achieved within the jth lower-level eqBound entry
  subPhiPREQCounts{k} = subPhiBREQCounts{k}(indPREQ2subBREQ(:,k),:);
end

% Series of intermediate bookkeeping steps to compute all possible
% combinations of energies between compatible partitioned groups
W1 = subPhiPREQ{1}';
W2 = subPhiPREQ{2}';
C1 = subPhiPREQCounts{1}';
C2 = subPhiPREQCounts{2}';

[nW,mW] = size(W1);

if nW*mW < maxSize % too small, don't divide
  
  C12 = kron(C1,ones(size(C2,1),1)).*kron(ones(size(C1,1),1),C2);
  W1 = full(W1); W2 = full(W2);
  W1(W1==0) = nan; W1 = W1 - subPhiBREQOffset(1);
  W2(W2==0) = nan; W2 = W2 - subPhiBREQOffset(2);
  WK1 = kron(W1,ones(size(W2,1),1));
  WK2 = kron(ones(size(W1,1),1),W2);
  WK12 = WK1 + WK2 + repmat(phiInterGroup',size(WK1,1),1);
  
  % uniqPhiPREQ = vector of all unique feasible energy values at current
  % level of subdivision
  [uniqPhiPREQ,~,uniqPhiPREQInds] = unique(WK12(~isnan(WK12)));
  uniqPhiPREQ = reshape(uniqPhiPREQ,length(uniqPhiPREQ),1);
  [~,WQ] = find(~isnan(WK12));
  
  % phiPREQCounts(i,j) = number of ways the unique energy value
  % uniqPhiPREQ(i) is achieved within the jth eqPart entry
  try 
    phiPREQCounts = sparse(uniqPhiPREQInds,WQ,C12(C12>0),length(uniqPhiPREQ),size(WK12,2))';
  catch e
    rethrow(e);
  end
else
  
  nDiv = ceil(nW*mW / maxSize);
  divInds = round(linspace(1,mW+1,nDiv+1));
  
  % Divide into pieces for parallel processing
  W1div = cell(nDiv,1); W2div = cell(nDiv,1);
  C1div = cell(nDiv,1); C2div = cell(nDiv,1);
  phiInterGroupDiv = cell(nDiv,1);
  for i = 1:nDiv
    W1div{i} = full(W1(:,divInds(i):divInds(i+1)-1));
    W2div{i} = full(W2(:,divInds(i):divInds(i+1)-1));
    W1div{i}(W1div{i}==0) = nan; W1div{i} = W1div{i} - subPhiBREQOffset(1);
    W2div{i}(W2div{i}==0) = nan; W2div{i} = W2div{i} - subPhiBREQOffset(2);
    C1div{i} = C1(:,divInds(i):divInds(i+1)-1);
    C2div{i} = C2(:,divInds(i):divInds(i+1)-1);
    phiInterGroupDiv{i} = phiInterGroup(divInds(i):divInds(i+1)-1)';
  end
  
  % Compute all feasible energy sums in parallel
  uniqPhiPREQDiv = cell(nDiv,1);
  Idiv = cell(nDiv,1);
  Jdiv = cell(nDiv,1);
  Vdiv = cell(nDiv,1);
  parfor i = 1:nDiv
    C12 = kron(C1div{i},ones(size(C2div{i},1),1)).*kron(ones(size(C1div{i},1),1),C2div{i});
    WK1 = kron(W1div{i},ones(size(W2div{i},1),1));
    WK2 = kron(ones(size(W1div{i},1),1),W2div{i});
    WK12 = WK1 + WK2 + repmat(phiInterGroupDiv{i},size(WK1,1),1);
    
    [uniqPhiPREQDiv{i},~,uniqPhiPREQInds] = unique(WK12(~isnan(WK12)));
    [~,WQ] = find(~isnan(WK12));
    
    Idiv{i} = uniqPhiPREQInds;
    Jdiv{i} = WQ + divInds(i) - 1;
    Vdiv{i} = C12(C12>0);
  end
  
  uniqPhiPREQ = unique(cat(1,uniqPhiPREQDiv{:}));
  for i = 1:nDiv    
    uniqPhiDiv2Full = find(ismember(uniqPhiPREQ,uniqPhiPREQDiv{i}));
    Idiv{i} = uniqPhiDiv2Full(Idiv{i});   
  end
  
  I = cat(1,Idiv{:});
  J = cat(1,Jdiv{:});
  V = cat(1,Vdiv{:});
  
  phiPREQCounts = sparse(I,J,V,length(uniqPhiPREQ),mW)';
end

