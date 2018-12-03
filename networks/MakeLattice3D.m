function [Adj,Vxyz] = MakeLattice3D(nrows,ncols,nslices)

n = nrows*ncols*nslices;
Adj = zeros(n,n);
nodesPerSlice = nrows*ncols;
Z = (1/nslices:1/nslices:1) - 1/(2*nslices);

for k = 1:nslices
  startInd = nodesPerSlice*(k-1) + 1;
  endInd = nodesPerSlice*k;
  [Ak,Vxyk] = MakeLattice2D(nrows,ncols);
  Adj(startInd:endInd,startInd:endInd) = Ak;
  Vxyz(startInd:endInd,:) = [Vxyk Z(k)*ones(nodesPerSlice,1)];
  
  if k < nslices
    Adj(startInd:endInd,startInd + nodesPerSlice:endInd + nodesPerSlice) = eye(nodesPerSlice);
    Adj(startInd + nodesPerSlice:endInd+nodesPerSlice,startInd:endInd) = eye(nodesPerSlice);
  end
  if k > 1
    Adj(startInd:endInd,startInd - nodesPerSlice:endInd - nodesPerSlice) = eye(nodesPerSlice);
    Adj(startInd - nodesPerSlice:endInd - nodesPerSlice,startInd:endInd) = eye(nodesPerSlice);
  end
end



