function ndx = gridPartition(Vxy,k)
%
% function [ndx,Pi,cost]= grPartition(C,k,nrep);
%
% Partitions the n-node undirected graph G defined by the matrix C
% 
% Inputs:
% A - n by n edge-weights matrix. In particular, c(i,j)=c(j,i) is equal 
%     to the cost associated with cuting the edge between nodes i and j.
%     This matrix should be symmetric and doubly stochastic. If this
%     is not the case, this matrix will be normalized to
%     satisfy these properties (with a warning).
% 
% Outputs:
% ndx  - n-vector with the cluster index for every node 
%       (indices from 1 to k)

xDiv = floor(sqrt(k));
yDiv = floor(k/xDiv);

X = Vxy(:,1);
Y = Vxy(:,2);

xMax = max(X);
xMin = min(X);
xRange = xMax-xMin;
yMax = max(Y);
yMin = min(Y);
yRange = yMax-yMin;

xDivs = xMin:xRange/xDiv:xMax;
xDivs(end) = inf;
yDivs = yMin:yRange/yDiv:yMax;
yDivs(end) = inf;

n = length(X);
ndx = zeros(n,1);

group = 1;
for i = 1:xDiv,
  xiMin = xDivs(i);
  xiMax = xDivs(i+1);
  xInd = X >= xiMin & X < xiMax;
  for j = 1:yDiv,
    yiMin = yDivs(j);
    yiMax = yDivs(j+1);
    groupInd = xInd & Y >= yiMin & Y < yiMax;
    ndx(groupInd) = group;
    group = group + 1;
  end
end

end

    
    
        




