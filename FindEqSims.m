% Find all equilibria for a range of network types / parameters

addpath(genpath('./'));

% Simulation parameters
maxGroupSize = 4; % q
numGroups = 2; % p
seed = 3;
numReps = 1;
basicLimit = 60;
exLimit = 15;

computeEnergies = true;
runBasic = true;
parallel = false;
saveResults = true;
drawThings = false;

networkType = 'tree';

switch networkType
  case 'torus'
    sizes = [(10:1:20) (25:5:65) (70:10:100) 250 500 1000];
    runExhaustive = true;
  case 'tree'
    sizes = (1:11);
    runExhaustive = false;
  case 'geometric'
    sizes = round(2.^(3:.5:5));
    runExhaustive = false;
end

numSizes = length(sizes);

if parallel
  pool = gcp('nocreate');
  if isempty(pool)
    parpool('local');
  end
end

meanLocalTime = zeros(numSizes,1);
meanGroupTime = zeros(numSizes,1);
meanExpandTime = zeros(numSizes,1);
meanTotalTime = zeros(numSizes,1);
meanCompBasic = zeros(numSizes,1);
meanCompExhaustive = zeros(numSizes,1);
meanNumEqs = zeros(numSizes,1);
meanNumOptEqs = zeros(numSizes,1);

for run = 1:numSizes
  
  localTime = zeros(numReps,1);
  groupTime = zeros(numReps,1);
  expandTime = zeros(numReps,1);
  totalTime = zeros(numReps,1);
  compBasic = zeros(numReps,1);
  compExhaustive = zeros(numReps,1);
  numEqs = zeros(numReps,1);
  numOptEqs = zeros(numReps,1);
 
  for rep = 1
    
    switch networkType
      case 'geometric'
        expDegree = 5;
        connectRange = sqrt((1+expDegree)/(pi*networkSize));
        [Adj,Vxy,Edges] = CreateGraph('type',networkType,'n',networkSize,...
          'connectRange',connectRange,'connected',true,'seed',rep);
      case 'torus'
        networkSize = sizes(run);
        nrows = sizes(run);
        ncols = 1;
        if strcmp(networkType,'torus')
          networkSize = ncols*nrows;
        end
        [Adj,Vxy,Edges] = CreateGraph('type',networkType,'nrows',nrows,...
          'ncols',ncols);
      case 'tree'
        childProbs = [0 0 1];
        numGen = sizes(run);
        [Adj,~,Edges] = CreateGraph('type','tree','childProbs',childProbs,...
          'numGen',numGen,'seed',rep + seed);
        if drawThings
          G = graph(Adj); fig = figure; GP = plot(G,'layout','force');
          Vxy = [GP.XData; GP.YData]'; close(fig);
        else
          n = size(Adj,1);
          Vxy = rand(n,2);
        end
        sizes(run) = n;
    end
    n = size(Adj,1); m = size(Edges,1);
    
    % For now, we assign all agents the same payoff matrix
    payoffMatrix = [0 3;1 0];
    nodePayoffs = repmat(payoffMatrix,[1 1 n]);
    
    % We also assume all agents use the best-response update rule
    updateRules = repmat({'FastBR'},n,1);
    
    % Energy/cost function
    energyFunction = @IsingEnergy;
    
    % Construct an asynchronous network game from this network and game dynamics
    ng = netgame('Adj',Adj,'Vxy',Vxy,'Edges',Edges,'updateRules',updateRules,...
      'payoffs',nodePayoffs,'networkSize',n,'synchronous',false,...
      'seed',0,'payoffSource','nodes','showGames',false,'initType',...
      'random','pauseGames',false,'saveFigures',false,...
      'networkType',networkType,'energyFunction',energyFunction);
    
    % Compute local equilibria
    tic;
    FindLocalEq(ng);
    localTime(rep) = toc;
    
    if runBasic && n <= basicLimit
      tic;
      JigsawGroups(ng,numGroups,inf,(1:ng.n)',computeEnergies,drawThings);
      basicTime = toc;
      compBasic(rep) = localTime(rep) + basicTime;
    end
    
    if runExhaustive && n <= exLimit
      tic;
      allEq = FindEqExhaustive(ng);
      exTime = toc;
      compExhaustive(rep) = exTime;
    end
    
    % Compute group equilibria
    tic;
    GroupData = JigsawGroups(ng,numGroups,maxGroupSize,(1:ng.n)',computeEnergies,drawThings);
    groupTime(rep) = toc;
    
    numEqs(rep) = GroupData.numEq;
    
    if computeEnergies
      maxPhi = max(GroupData.uniqPhiPREQ);
      % Compute optimal equilibria
      [eqOpt,compOpt] = ExpandEq(GroupData);
      expandTime(rep) = toc;
      numOptEqs(rep) = size(eqOpt,1);
    end
    
    totalTime(rep) = localTime(rep) + groupTime(rep) + expandTime(rep);
    
  end
  
  meanLocalTime(run) = mean(localTime);
  meanGroupTime(run) = mean(groupTime);
  meanExpandTime(run) = mean(expandTime);
  meanTotalTime(run) = mean(totalTime);
  meanCompBasic(run) = mean(compBasic);
  meanCompExhaustive(run) = mean(compExhaustive);
  meanNumEqs(run) = mean(numEqs);
  meanNumOptEqs(run) = mean(numOptEqs);
  
  figure(2); clf;
  semilogy(sizes,meanCompExhaustive,'s-','Color',[.7 .2 .2],'LineWidth',2); hold on;
  semilogy(sizes,meanCompBasic,'d-','Color',[.2 .7 .2],'LineWidth',2);
  semilogy(sizes,meanTotalTime,'o-','Color',[.2 .2 .7],'LineWidth',2); hold on;
  set(gca,'FontSize',16);
  xlabel('network size');
  ylabel('computation time (s)');
  legend('Exhaustive Search','Local Intersections','Recursive Partitioning','Location','SouthEast');
  grid on; box on;
  drawnow;
  
end




