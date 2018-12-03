function L = LogMat(n,allS)
% Generate logical matrix containing all possible configurations of n
% strategy states as rows of the matrix
ns = length(allS);

L = zeros(ns^n,n);
for i = 1:n
  q = allS;
  for j = 1:n-i
    q = kron(q,ones(size(allS)));
  end
  L(:,i) = repmat(q,ns^(i-1),1);
end
