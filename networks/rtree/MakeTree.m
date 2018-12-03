function [A,Vxy,parents,genNodes]=MakeTree(ngen, prob_vec, seed, minNodes, maxNodes)
% BRANCH Simulate a branching process and plot it as a tree.
%   Each individual can get 0-N children with probabilities
%   p0-pN.
%
% [parents] = branch(ngen, prob_vec)
%
% Inputs:
%   ngen - number of generations to proceed
%   prob_vec - vector of probabilities for children. Note
%     that the first entry is a probability for 0.
%
% Outputs:
%   parents - a random tree in a random order. parents(i) is the index
%     number of the parents to the node i.
%

% Authors: R.Gaigalas, I.Kaj
% v1.5 Created 07-Nov-01
%      Modified 23-Nov-05 changed variable names

if (nargin<1) % default parameter values
  ngen = 6;
  prob_vec = [0.2 0.3 0.1 0.4];
end
if nargin < 4,
  minNodes = 0;
end
if nargin < 5
  maxNodes = inf;
end

bad = true;

while bad, % Loop until tree meets specifications  
  rng(seed);
  
  max_kids = length(prob_vec)-1; % maximal number of children
  cum_prob = [cumsum(prob_vec) 1]; % distribution function
  parents = 0;
  extinct = 0;
  nkids = 1; % number of children in the last generation
  genNodes = cell(ngen+1,1);
  genNodes{1} = 1;
  
  for i=1:ngen
    
    % indices of individuals in the previous generation
    gprevi = length(parents)-nkids+1:length(parents);
    npar = nkids;
    
    % keep track of nodes by generation
    if i>1,
      genNodes{i} = gprevi;
    end
    
    % generate a set of numbers of children for all individuals
    % according to the distribution
    runi = rand(1, npar);  % a uniform sample
    
    % add the kids to the tree
    nkids = 0;
    for j=1:max_kids
      % find indexes of the parents who give birth to j children
      jpi = gprevi((runi>cum_prob(j)) ...
        & (runi<=cum_prob(j+1)));
      if (~isempty(jpi))
        parents = [parents repmat(jpi, 1, j)];
        nkids = nkids+length(jpi)*j;
      end
    end
    
    % exit if no individuals left
    if (nkids==0)
      extinct = 1;
      break;
    end
  end
  
  % Add last generation nodes to data cell
  numParentNodes = 0;
  for i = 1:ngen,
    numParentNodes = numParentNodes + length(genNodes{i});
  end
  n = length(parents);
  genNodes{ngen+1} = numParentNodes+1:n;
  
%   if (extinct)
%     fprintf('Extinct in generation %d\n', i);
%   else
%     fprintf('Stopped in generation %d\n', ngen);
%   end
  
  % figure(1);
  [x,y,~]=trimtreelayout(parents);
  Vxy = [x' y'];
  % figure(2);
  % newtreeplot(parents);
  
  if length(parents) >= minNodes && length(parents) <= maxNodes,
    bad = false;
  end
  
  seed = seed + 10000;
end

A = zeros(n,n);
for i = 2:length(parents),
  A(i,parents(i)) = 1;
  A(parents(i),i) = 1;
end

