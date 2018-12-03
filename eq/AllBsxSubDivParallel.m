function out = AllBsxSubDivParallel(func,A,B)
% Subdivide a call to bsxfun if the generated array size will be too big

% Largest size
maxSize = 1e9;

% Compute size of largest matrix
nA = size(A,3);
nBox = size(B,1)*numel(A);

if nBox < maxSize % too small, don't bother
  Z = all(bsxfun(func,A,B),2);
  out = sparse(reshape(Z,size(B,1),size(A,3)));
else
  
  nDiv = ceil(nBox / maxSize);
  divInds = round(linspace(1,nA+1,nDiv+1));
  
  % Divide into pieces for parallel processing
  Aslice = cell(nDiv,1);
  for i = 1:nDiv
    Aslice{i} =  A(:,:,divInds(i):divInds(i+1)-1);
  end
  
  outSlice = cell(nDiv,1);
  parfor i = 1:nDiv
    display(['AllBsxSubDiv: ' num2str(i) ' of ' num2str(nDiv)]);
    Z = all(bsxfun(func,Aslice{i},B),2);
    outSlice{i} = sparse(reshape(Z,size(B,1),size(Aslice{i},3)));
  end
  
  out = cat(2,outSlice{:});
  
end