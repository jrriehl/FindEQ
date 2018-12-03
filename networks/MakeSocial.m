function [Adj,Vxy,Edges] = MakeSocial(n,seed)

r0 = 0.005; % nominal meeting rate
r1 = .1; % preferential meeting rate
gamma = 0.005; % decay rate
Adj = zeros(n,n);
Edges = zeros(0,2);
maxDeg = 4;
deg = zeros(n,1);
numIter = 100;
fracMax = 0.95;

rng(seed);
allPairs = [kron((1:n)',ones(n,1)) kron(ones(n,1),(1:n)')];
allPairs(allPairs(:,1)>=allPairs(:,2),:) = [];
np = size(allPairs,1);

init = false;

for j = 1:numIter
  
  % 1. Choose random pairs
  chooseInds = randperm(np,round(np*r0));
  choosePairs = allPairs(chooseInds,:);
  
  % Add any of these that don't already exist and don't violate maximum degree
  duplicates = ismember(choosePairs,Edges,'rows');
  maxDegNodes = find(deg >= maxDeg);
  keepInds = ~duplicates & ~any(ismember(choosePairs,maxDegNodes),2);
  keepPairs = choosePairs(keepInds,:);
  
  Edges = [Edges;keepPairs];
  Adj = full(sparse(Edges(:,1),Edges(:,2),ones(size(Edges,1),1),n,n));
  Adj = (Adj + Adj');
  deg = sum(Adj,2);
  
  % 2. Choose random vertices with probabilities proportional to deg(deg-1)
  numPairNeighbors = (deg.*(deg-1))/2;
  numChooseV = round(sum(numPairNeighbors)*r1);
  for i = 1:numChooseV
    posV = find(numPairNeighbors);
    if isempty(posV), continue; end
    probs = numPairNeighbors(posV)/sum(numPairNeighbors(posV));
    cumProbs = cumsum(probs);
    r = rand;
    selInd = find(cumProbs > r,1,'first');
    selNode = posV(selInd);
    
    % Choose random pair of those nodes neighbors to connect
    neighbors = find(Adj(:,selNode));
    nn = length(neighbors);
    allNeighborPairs = [kron(neighbors,ones(nn,1)) kron(ones(nn,1),neighbors)];
    allNeighborPairs(allNeighborPairs(:,1)>=allNeighborPairs(:,2),:) = [];
    nnp = size(allNeighborPairs,1);
    
    selPair = allNeighborPairs(ceil(rand*nnp),:);
    maxDegNodes = find(deg >= maxDeg);
    if ~ismember(selPair,Edges,'rows') && ~any(ismember(selPair,maxDegNodes))
      Edges = [Edges;selPair];
    end
    
    numPairNeighbors(selNode) = 0;
    Adj = full(sparse(Edges(:,1),Edges(:,2),ones(size(Edges,1),1),n,n));
    Adj = (Adj + Adj');
    deg = sum(Adj,2);
  end
  
  if mean(deg >= 5) >= fracMax
    init  = false;
  end
  
  % 3. Randomy decay some connections and delete
  if ~init
    numEdges = size(Edges,1);
    numChooseV = min(n,ceil(rand*gamma*numEdges));
    chooseV = randperm(n,numChooseV);
    for i = 1:numChooseV
      probs = deg(chooseV);
      cumProbs = cumsum(probs);
      r = rand;
      selInd = find(cumProbs > r,1,'first');
      selNode = chooseV(selInd);
      
      neighbors = find(Adj(:,selNode));
      if isempty(neighbors), continue; end
      selNeighborInd = ceil(rand*length(neighbors));
      selNeighbor = neighbors(selNeighborInd);
      Edges(ismember(Edges,[selNode selNeighbor],'rows'),:) = [];
      Adj = full(sparse(Edges(:,1),Edges(:,2),ones(size(Edges,1),1),n,n));
      Adj = (Adj + Adj');
      deg = sum(Adj,2);
    end
  end
  
end

G = graph(Adj);
GP = plot(G,'layout','force');
Vxy = [GP.XData; GP.YData]';
[I,J] = find(Adj);
Edges = [I J];



