function [A,Vxy,E,parents,genNodes] = CreateGraph(varargin)
%CREATEGRAPH - constructs several different types of graphs and returns the
%adjecency matrix A and vertex positions Vxy for plotting
%
% Output arguments parents and genNodes only used for trees

% Default values
minNodes = 0;
maxNodes = inf;
parents = 0;
genNodes = 0;
forceConnect = false;

if (rem(nargin,2) == 1)
  error('Number of arguments should be even') ;
else
  for i=1:2:nargin
    arg1 = varargin{i};
    arg2 = varargin{i+1};
    switch(arg1)
      case 'type'
        type = arg2;
      case 'n'
        n = arg2;
      case 'nrows'
        nrows = arg2;
      case 'ncols'
        ncols = arg2;
      case 'numGen'
        numGen = arg2;
      case 'childProbs'
        childProbs = arg2;
      case 'maxNodes'
        maxNodes = arg2;
      case 'minNodes'
        minNodes = arg2;
      case 'rewireProb'
        rewireProb = arg2;
      case 'connectProb'
        connectProb = arg2;
      case 'connectRange'
        connectRange = arg2;
      case 'neighborDist'
        neighborDist = arg2;
      case 'minDegree'
        minDegree = arg2;
      case 'connected'
        forceConnect = arg2;
      case 'seed'
        seed = arg2;
        rng(seed);
    end
  end
end

disconnected = true;
% If forceConnect = true, repeat until the randomly generated graph is connected
while disconnected
  
  switch(type)
    case 'lattice'
      [A,Vxy] = MakeLattice2D(nrows,ncols);
    case 'torus'
      [A,Vxy] = MakeTorus(nrows,ncols);
    case 'star'
      childProbs = zeros(1,n); childProbs(end) = 1;
      [A,Vxy] = MakeTree(1,childProbs,0,0,inf);
    case 'tree'
      [A,Vxy,parents,genNodes] = MakeTree(numGen,childProbs,seed,minNodes,maxNodes);      
    case 'erdos-renyi'
      [A,Vxy] = MakeRandom(n,connectProb);
    case 'small-world'
      A = smallw(n,neighborDist,rewireProb);
    case 'scale-free' %preferential attachment
      A = pref(n,minDegree);
    case 'geometric'
      [A,Vxy] = MakeRandGeo(n,connectRange);
    case 'dir-geometric'
      [A,Vxy] = MakeRandDirGeo(n,connectRange);
    case 'complete'
      A = ones(n) - eye(n);
  end
  
  if forceConnect
    % Check connectivity
    D = diag(sum(A));
    L = D - A;
    lambda = eig(L);
    if lambda(2) < 1e-8
      disconnected = true;
    else
      disconnected = false;
    end
  else
    disconnected = false;
  end
  
end

if ~exist('Vxy','var') || strcmp(type,'torus') && (nrows==1 || ncols==1)
  % Put vertices on a circle
  n = size(A,1);
  theta = 2*pi*(1:n)'/n;
  Vxy = [cos(theta) sin(theta)];
end

% Edges (directed)
[I,J] = find(A); E = [J I];

end