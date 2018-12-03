function [parents]=branch(ngen, prob_vec)
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
 
 max_kids = length(prob_vec)-1 % maximal number of children
 cum_prob = [cumsum(prob_vec) 1]; % distribution function
 parents = 0; 
 extinct = 0;
 nkids = 1; % number of children in the last generation
 
 for i=1:ngen
   
   % indices of individuals in the previous generation
   gprevi = length(parents)-nkids+1:length(parents);
   npar = nkids;

   % generate a set of numbers of children for all individuals
   % according to the distribution
   runi = rand(1, npar);  % a uniform sample
   
   % add the kids to the tree   
   nkids = 0;
   for j=1:max_kids
     % find indexes of the parents who give birth to j children 
     jpi = gprevi(find((runi>cum_prob(j)) ...
                     & (runi<=cum_prob(j+1))));
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
 
 if (extinct)
   fprintf('Extinct in generation %d\n', i);
 else
   fprintf('Stopped in generation %d\n', ngen);
 end
 
 figure(1);
 Vxy = trimtreeplot(parents);   
% figure(2);
% newtreeplot(parents);    

n = length(parents);
A = zeros(n,n);
for i = 1:length(parents),
  A(i,parents(i)+1) = 1;
  A(parents(i)+1,i) = 1;
end
 
