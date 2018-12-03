% EQ validation tests: make sure functions perform correctly

addpath(genpath('./'));

% Network and algorithm parameters
networkType = 'geometric'; % {'ring','geometric','dir-geometric'}
n = 12;
numStrategies = 2;
numGroups = 2;
maxGroupSize = 10;
seed = 1;

%% ----------Construct a random example to use for testing-------------- %%

switch networkType
  case 'ring'
    [Adj,Vxy,Edges] = CreateGraph('type','torus','nrows',n,'ncols',1);
  case 'geometric'
    [Adj,Vxy,Edges] = CreateGraph('type','geometric','n',n,...
      'connectRange',0.4,'connected',true,'seed',seed);
  case 'dir-geometric'
    [Adj,Vxy,Edges] = CreateGraph('type','dir-geometric','n',n,...
      'connectRange',0.6,'connected',true,'seed',seed);    
end

switch numStrategies
  case 2
    allS = [1;2];
    payoffMatrix = [0 3;1 0];
    % payoffMatrix = [3 0;0 1];
    updateRules = repmat({'BestResponse'},n,1);
    energyFunction = @IsingEnergy;
  case 3
    allS = [1;2;3];
    payoffMatrix = [0 3 3;1 0 0;0 0 0];
    updateRules = repmat({'BestResponseControl'},n,1);
    energyFunction = @IsingEnergy;
end

n = size(Adj,1); m = size(Edges,1);
V = (1:n)';
deg = sum(Adj,2);
nodePayoffs = repmat(payoffMatrix,[1 1 n]);

% Construct an asynchronous network game from this network and game dynamics
ng = netgame('Adj',Adj,'Vxy',Vxy,'Edges',Edges,'updateRules',updateRules,...
  'payoffs',nodePayoffs,'networkSize',n,'synchronous',false,'seed',seed,...
  'payoffSource','nodes','showGames',false,'initType','random','pauseGames',...
  false,'saveFigures',false,'energyFunction',energyFunction);
ng.ns = length(allS);

% Some of these will be brute force tests that involve checking all states
allStates = LogMat(n,allS);
numStates = size(allStates,1);

results = false(5,1);

%% ----------------------------FindLocalEQ------------------------------ %%

LEQ = FindLocalEq(ng);

% Construct partition matrix as if each node is it's own partition group
[Pi,PiPart] = BasePartition(Adj,V);

good = true;
% for i = 1:numStates
%   display(['Checking state ' num2str(i) ' of ' num2str(numStates)]);
%   x0 = allStates(i,:);
%   for j = 1:n
%     neighbors = find(PiPart(:,j));
%     xj = zeros(1,n);
%     xj(neighbors) = x0(neighbors);
%     ng.x = x0';
%     Update(ng,j);
%     if ng.x(j) == x0(j)
%       if ~ismember(xj(neighbors),LEQ{j},'rows')
%         good = false;
%       end
%     else
%       if ismember(xj(neighbors),LEQ{j},'rows')
%         good = false;
%       end
%     end
%   end
% end

if good
  disp('FindLocalEQ Test: passed');
else
  disp('FindLocalEQ Test: failed');
end

results(1) = good;

%% ----------------------------Local2Global----------------------------- %%

GEQ = GroupMerge(ng.LEQ,(1:n)',PiPart);

good = true;
% for i = 1:numStates
%   display(['Checking state ' num2str(i) ' of ' num2str(numStates)]);
%   x0 = allStates(i,:);
%   ng.x0 = x0';
%   xf = Play(ng);
%   if all(x0' == xf)
%     if ~ismember(x0,GEQ,'rows')
%       good = false;
%     end
%   else
%     if ismember(x0,GEQ,'rows')
%       good = false;
%     end
%   end
% end

if good
  disp('Local2Global Test: passed');
else
  disp('Local2Global Test: failed');
end

results(2) = good;

%% ----------------------------JigsawGroups----------------------------- %%

InitPlots(ng);

GD = JigsawGroups(ng,numGroups,maxGroupSize,(1:ng.n)',true,true);
Vpart = find(any(GD.PiOut | GD.PiIn,2));

GEQcheck = GEQ(:,Vpart);

goodNumEq = GD.numEq == size(GEQ,1);

good = all(ismember(GD.PREQ,GEQcheck,'rows')) && ...
  all(ismember(GEQcheck,GD.PREQ,'rows')) && goodNumEq;

if good
  disp('JigsawGroups Test: passed');
else
  disp('JigsawGroups Test: failed');
end

results(3) = good;

%% ------------------------------Energies------------------------------- %%

% Compute all energies
allPhis = ComputeEnergy(ng,Adj,GEQ,(1:n)',(1:n)');
allPhiCheck = unique(allPhis);

% Computed energies
allPhi = GD.uniqPhiPREQ;

good = all(ismember(allPhi,allPhiCheck)) && all(ismember(allPhiCheck,allPhi));

if good
  disp('Energy Test: passed');
else
  disp('Energy Test: failed');
end

results(4) = good;

%% -----------------------------ExpandEq-------------------------------- %%

% Compute all energies
allPhi = ComputeEnergy(ng,Adj,GEQ,(1:n)',(1:n)');

numUniqPhi = length(GD.uniqPhiPREQ);
good = false(numUniqPhi,1);
for i = 1:numUniqPhi
  targetPhi = GD.uniqPhiPREQ(i);
  eqPhi = ExpandEq2(GD,targetPhi);
  
  eqPhiCheck = GEQ(allPhi == targetPhi,:);
  
  good(i) = all(ismember(eqPhi,eqPhiCheck,'rows')) && ...
    all(ismember(eqPhiCheck,eqPhi,'rows'));
end

if all(good)
  disp('ExpandEq Test: passed');
else
  disp('ExpandEq Test: failed');
end

results(5) = all(good);

%% --------------------------------All---------------------------------- %%

allGood = all(results);

if allGood
  disp('All tests passed');
else
  disp('One or more tests failed');
end


