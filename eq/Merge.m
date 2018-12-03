function X = Merge(X1,X2,V1,V2)
% Compute the row-intersection of two matrices X1 and X2, where zero
% entries are assumed to be unknown (can take any value). This function
% uses full matrices (to save memory, believe it or not).
%
% For example, in a binary system (x in {1,2}), the row X1 = [0 1 2]  
% represents the set X1 = {[1 1 2],[2 1 2]}. And the row X2 = [1 0 2]
% represents the set X2 = {[1 1 2],[1 2 2]}. The intersection of X1 and
% X2 is thus X = {[1 1 2]}. Zero entries that are common to both X1 and X2
% are preserved in the about because there is no new information. If the
% intersection of X1 and X2 is empty, then X needs to represent all
% possible combinations of the entries of X1 and X2.

Vintersect = V1 & V2;
Vunion = V1 | V2;
Vdiff = V2 & ~V1;
V1ind = find(V1);
V2ind = find(V2);
VIind = find(Vintersect);
VUind = find(Vunion);
VDind = find(Vdiff);
X1a = X1(:,ismember(V1ind,VIind));
try
  X2a = X2(:,ismember(V2ind,VIind));
catch e
  rethrow(e);
end

% TODO: Fix this part
% In the (hopefully rare) case of no overlapping columns, we just need to 
% make a stack with all permutations of states between X1 and X2
if isempty(VIind)
  stack = @(x,y) [kron(x,ones(size(y,1),1)) kron(ones(size(x,1),1),y)];
  V12 = [find(V1);find(V2)];
  [~,V12ord] = sort(V12);
  X1v = double(X1);
  X2v = double(X2);
  nrows = size(X1v,1)*size(X2v,1);
  ncols = sum(V1) + sum(V2);
  I = repmat((1:nrows)',1,ncols);
  J = repmat((1:ncols),nrows,1);
  V = stack(X1v,X2v);
  X = uint8(full(sparse(I,J,V,nrows,ncols)));
  X = X(:,V12ord);
  return;
end

% Otherwise (normally), check for matches between X1 and X2
matches = AllBsxSubDivParallel(@eq,reshape(X2a.',1,length(VIind),[]),X1a);
[I,J] = find(matches); I = reshape(I,length(I),1); J = reshape(J,length(J),1);
X = zeros(size(I,1),length(VUind),'uint8');
X(:,ismember(VUind,V1ind)) = X1(I,:);
X(:,ismember(VUind,VDind)) = X2(J,ismember(V2ind,VDind));