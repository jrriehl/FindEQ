function out = IsingEnergy(Adj,X,Xind,selNodes)

X = X(:,ismember(Xind,selNodes));
AdjSel = Adj(selNodes,selNodes);

X = double(X);
X = 2*X - 3;

out = -sum((X*AdjSel.').*X,2);