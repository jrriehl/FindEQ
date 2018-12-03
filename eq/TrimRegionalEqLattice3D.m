function X = TrimRegionalEqLattice3D(ng,X,Vpart,Vout)

bad = false(size(X,1),1);

for ii = 1:length(Vout)
  i = find(ismember(Vpart,Vout(ii)));
  outNeighbors = find(ng.Adj(Vout(ii),:));
  neighborInds = ismember(Vpart,outNeighbors);
  neighborStates = X(:,neighborInds);
  
  one = X(:,i) == 1;
  two = X(:,i) == 2;
  
  numOneNeighbors = sum(neighborStates == 1,2);
  numTwoNeighbors = sum(neighborStates == 2,2);
  
  iBad = one & numOneNeighbors == 3 | two & numTwoNeighbors == 3;
  bad = bad | iBad;
end

X = X(~bad,:);