classdef netgame < handle
  % Class for modeling a network game, where the nodes correspond to
  % players and the edges correspond to 2-player matrix games between the
  % connected nodes. Allows arbitrary payoff matrices and various update
  % rules.
  properties
    Adj;              % graph adjacency matrix (n,n)
    Adjs;             % sparse adjacency matrix (n,n)
    Edges;            % 2-column edge list (m,2)
    Vxy;              % Planar vertex coordinates (n,2)
    payoffs;          % payoff matrices for all agents (ns,ns,n)
    edgePayoffs;      % payoff matrices for all edges (ns,ns,m)
    payoffSource;     % where payoffs are defined {'nodes','edges'}
    updateRule;       % network update method
    updateRules;      % update rule for each agent
    synchronous;      % 0=asynchronous, 1=synchronous, 1/2=partial sync
    n;                % number af agents in network
    m;                % number of edges in network
    ns;               % number of available strategies
    x0;               % initial strategy state (n,1)
    xd;               % desired strategy state (n,1)
    x;                % current strategy state (n,1)
    u;                % strategy control input (n,1)
    y;                % current agent payoffs  (n,1)
    h;                % characteristic payoff quantity (n,1)
    phi;              % local energy function (n,1)
    deltaA;           % payoff advantage to use A against A on a given edge (m,1)
    deltaB;           % payoff advantage to use B against B on a given edge (m,1)
    edgeAdjMap;       % adjacency matrix containing edge indices
    weightedAdj;      % adjacency matrix containing euclidian edge lengths
    connectPairs;     % set of all possible neighbor-pairs (only for TSP-problems)
    Pi;               % Stacked payoff matrix data (n,n,ns)
    controlNodes;     % agents whose strategy in controlled
    rewardBudget;     % budget limit for incentives
    terminated;       % true if game reaches equilibrium
    numRounds;        % max. number of game rounds before stopping
    lambda;           % proportional imitation rate
    logitBeta;        % exponential factor for logit response function
    deg;              % degree vector (n,1)
    LEQ;              % set of local equilibria at each node {n,1}
    LEQSizes;         % size of LEQ sets at each node (n,1)
    networkType;      % type of network {complete,ring,torus,geometric, etc.}
    networkSize;      % target size for creating network
    payoffType;       % type of payoff matrices
    coordLevel;       % base level of coordination for agents (default=1)
    coordFraction;    % fraction of coordinators when payoffs are split
    payoffVariance;   % amount of random variation in agent payoffs
    initType;         % how to initialize strategies {unform,random,randomEq}
    energyFunction;   % cost or energy function to evaluate
    connectRange;     % connection radius for geo. rand. networks
    neighborDist;     % number of neighbors to connect on each side for small-world
    rewireProb;       % small-world rewire probability
    minDegree;        % scale-free minimum degree parameter
    nrows;            % number of rows in lattice/torus network
    ncols;            % number of columns in lattice/torus network
    seed;             % random number generator seed
    hAgents;          % graphics handle for agents
    hNormal;          % graphics handle for uncontrolled agents
    hActive;          % graphics handle for active agents
    hSwitch;          % graphics handle for switching agents
    hText;            % graphics handle for agent labels
    nodeSize;         % node size for graph plot
    fast;             % if true, only activate agents will change
    showGames;        % if true, plot progress of games
    showStatus;       % if true, display progress of algorithm at various stages
    showActivations;  % if true, show active and switching agents
    pauseGames;       % if true, pause between each update round
    showEnergy;       % if true, display energy values with each update
    nodeLabels;       % true to display labels, false to display no text
    saveFigures;      % if true, save snapshots of each stage of the game
    avgPayoffs;       % if true, compute average payoffs instead of cumulative
    saveStr;          % root name of files to save
  end
  methods
    function ng = netgame(varargin)
      
      % Default game setup
      ng.networkSize = 15;
      ng.Adj = [];
      ng.weightedAdj = [];
      ng.updateRule = 'BestResponse';
      ng.updateRules = [];
      ng.synchronous = false;
      ng.terminated = false;
      ng.numRounds = 100;
      ng.networkType = 'none';
      ng.payoffType = 'coord';
      ng.payoffSource = 'nodes';
      ng.edgePayoffs = 0;
      ng.deltaA = 0;
      ng.deltaB = 0;
      ng.payoffVariance = 0;
      ng.coordLevel = 1;
      ng.coordFraction = 0.5;
      ng.ns = 2;
      ng.initType = 'uniform';
      ng.controlNodes = [];
      ng.rewardBudget = inf;
      ng.connectRange = 0.3;
      ng.neighborDist = 2;
      ng.rewireProb = 0.5;
      ng.nodeLabels = false;
      ng.nodeSize = 5;
      ng.minDegree = 2;
      ng.LEQ = [];
      ng.nrows = 5;
      ng.ncols = 5;
      ng.seed = 0;
      ng.lambda = 1;
      ng.logitBeta = 1;
      ng.fast = true;
      ng.connectPairs = [];
      ng.showGames = false;
      ng.showStatus = false;
      ng.showActivations = false;
      ng.pauseGames = true;
      ng.showEnergy = false;
      ng.saveFigures = false;
      ng.avgPayoffs = false;
      ng.saveStr = 'NetGame';
      
      %--------------------------------------------------------------------
      if (rem(nargin,2) == 1)
        error('Number of arguments should be even') ;
      else
        for i=1:2:nargin
          arg1 = varargin{i};
          switch(arg1)
            case 'Adj'
              ng.Adj = varargin{i+1};
              ng.n = size(ng.Adj,1);
            case 'Edges'
              ng.Edges = varargin{i+1};
              ng.m = size(ng.Edges,1);
            case 'Vxy'
              ng.Vxy = varargin{i+1};
            case 'payoffs'
              ng.payoffs = varargin{i+1};
              ng.ns = size(ng.payoffs,1);
              ng.payoffType = 'given';
            case 'payoffSource'
              ng.payoffSource = varargin{i+1};
            case 'edgePayoffs'
              ng.edgePayoffs = varargin{i+1};
              ng.ns = size(ng.edgePayoffs,1);
              ng.payoffType = 'given';
            case 'updateRule'
              ng.updateRule = varargin{i+1};
            case 'updateRules'
              ng.updateRules = varargin{i+1};
            case 'synchronous'
              ng.synchronous = varargin{i+1};
            case 'numRounds'
              ng.numRounds = varargin{i+1};
            case 'logitBeta'
              ng.logitBeta = varargin{i+1};
            case 'saveFigures'
              ng.saveFigures = varargin{i+1};
            case 'saveStr'
              ng.saveStr = varargin{i+1};
            case 'showGames'
              ng.showGames = varargin{i+1};
            case 'showActivations'
              ng.showActivations = varargin{i+1};
            case 'fast'
              ng.fast = varargin{i+1};
            case 'pauseGames'
              ng.pauseGames = varargin{i+1};
            case 'nodeLabels'
              ng.nodeLabels = varargin{i+1};
            case 'avgPayoffs'
              ng.avgPayoffs = varargin{i+1};
            case 'networkType'
              ng.networkType = varargin{i+1};
            case 'networkSize'
              ng.networkSize = varargin{i+1};
            case 'payoffType'
              ng.payoffType = varargin{i+1};
            case 'coordLevel'
              ng.coordLevel = varargin{i+1};
            case 'coordFraction'
              ng.coordFraction = varargin{i+1};
            case 'payoffVariance'
              ng.payoffVariance = varargin{i+1};
            case 'showEnergy'
              ng.showEnergy = varargin{i+1};
            case 'rewardBudget'
              ng.rewardBudget = varargin{i+1};
            case 'numStrategies'
              ng.ns = varargin{i+1};
            case 'initType'
              ng.initType = varargin{i+1};
            case 'nodeSize'
              ng.nodeSize = varargin{i+1};
            case 'energyFunction'
              ng.energyFunction = varargin{i+1};
            case 'neighborDist'
              ng.neighborDist = varargin{i+1};
            case 'rewireProb'
              ng.rewireProb = varargin{i+1};
            case 'minDegree'
              ng.minDegree = varargin{i+1};
            case 'nrows'
              ng.nrows = varargin{i+1};
            case 'ncols'
              ng.ncols = varargin{i+1};
            case 'seed'
              ng.seed = varargin{i+1};
          end
        end
      end
           
      
      % Create network
      if isempty(ng.Adj)
        switch ng.networkType
          case 'complete'
            [ng.Adj,ng.Vxy,ng.Edges] = CreateGraph('type','complete','n',ng.networkSize);
          case 'ring'
            [ng.Adj,ng.Vxy,ng.Edges] = CreateGraph('type','torus','nrows',ng.networkSize,'ncols',1);
          case 'torus'
            [ng.Adj,ng.Vxy,ng.Edges] = CreateGraph('type','torus','nrows',ng.nrows,'ncols',ng.ncols);
          case 'geometric'
            [ng.Adj,ng.Vxy,ng.Edges] = CreateGraph('type','geometric','n',ng.networkSize,...
              'connectRange',ng.connectRange,'connected',true,'seed',ng.seed);
          case 'dir-geometric'
            [ng.Adj,ng.Vxy,ng.Edges] = CreateGraph('type','dir-geometric','n',ng.networkSize,...
              'connectRange',ng.connectRange,'connected',true,'seed',ng.seed);
          case 'small-world'
            [ng.Adj,ng.Vxy,ng.Edges] = CreateGraph('type','small-world','n',ng.networkSize,...
              'neighborDist',ng.neighborDist,'rewireProb',ng.rewireProb,'seed',ng.seed);
          case 'scale-free'
            [ng.Adj,ng.Vxy,ng.Edges] = CreateGraph('type','scale-free','n',ng.networkSize,...
              'minDegree',ng.minDegree,'seed',ng.seed);
          case 'none'
            [ng.Adj,ng.Vxy,ng.Edges] = CreateGraph('type','complete','n',ng.networkSize);
        end
      end
      ng.Adjs = sparse(ng.Adj);
      ng.n = size(ng.Adj,1);
      ng.m = size(ng.Edges,1);
      
      % Set payoffs
      if strcmp(ng.payoffs ~= 0,'given')
        if strcmp(ng.payoffSource,'edges')
          % Set corresponding edge payoffs
          ng.edgePayoffs = zeros(ng.ns,ng.ns,ng.m);
          for e = 1:size(ng.Edges,1)
            i = ng.Edges(e,1);
            pi_ij = ng.payoffs(:,:,i);
            ng.edgePayoffs(:,:,e) = pi_ij;
          end
        end
      else
        if strcmp(ng.payoffSource,'nodes')
          switch ng.payoffType
            case 'coord'
              baseGame = eye(ng.ns);
              baseGame(1,1) = ng.coordLevel;
              ng.payoffs = repmat(baseGame,[1 1 ng.n]) + ng.payoffVariance*rand([ng.ns ng.ns ng.n]);
            case 'anticoord'
              baseGame = ones(ng.ns)-eye(ng.ns);
              baseGame(2,1) = ng.coordLevel;
              ng.payoffs = repmat(baseGame,[1 1 ng.n]) + ng.payoffVariance*rand([ng.ns ng.ns ng.n]);
            case 'random'
              ng.payoffs = rand([ng.ns ng.ns ng.n]);
            case 'split'
              ng.payoffs = zeros([ng.ns ng.ns ng.n]);
              numC = ceil(ng.n*ng.coordFraction);
              numAC = ng.n - numC;
              agentsC = randperm(ng.n,numC);
              agentsAC = setdiff(1:ng.n,agentsC);
              baseGameC = eye(ng.ns);
              baseGameC(1,1) = ng.coordLevel;
              baseGameAC = ones(ng.ns)-eye(ng.ns);
              baseGameAC(1,2) = ng.coordLevel;
              ng.payoffs(:,:,agentsC) = repmat(baseGameC,[1 1 numC]) + ng.payoffVariance*rand([ng.ns ng.ns numC]);
              ng.payoffs(:,:,agentsAC) = repmat(baseGameAC,[1 1 numAC]) + ng.payoffVariance*rand([ng.ns ng.ns numAC]);
            case 'row-coord'
              ng.payoffs = [2-cumsum(rand(1,2,ng.n));cumsum(rand(1,2,ng.n))];
          end
          ng.ns = size(ng.payoffs(:,:,1),1);
          
          % Set corresponding edge payoffs
          ng.edgePayoffs = zeros(ng.ns,ng.ns,ng.m);
          for e = 1:size(ng.Edges,1)
            i = ng.Edges(e,1);
            pi_ij = ng.payoffs(:,:,i);
            ng.edgePayoffs(:,:,e) = pi_ij;
          end
          
        else % define payoffs on edges
          switch ng.payoffType
            case 'coord'
              baseGame = eye(ng.ns);
              baseGame(1,1) = ng.coordLevel;
              ng.edgePayoffs = repmat(baseGame,[1 1 ng.m]) + ng.payoffVariance*rand([ng.ns ng.ns ng.m]);
            case 'anticoord'
              baseGame = ones(ng.ns) - eye(ng.ns);
              baseGame(1,:) = baseGame(1,:)*ng.coordLevel;
              ng.edgePayoffs = repmat(baseGame,[1 1 ng.m]) + ng.payoffVariance*rand([ng.ns ng.ns ng.m]);
            case 'random'
              ng.edgePayoffs = rand([ng.ns ng.ns ng.m]);
            case 'random-undirected'
              ng.edgePayoffs = zeros([ng.ns ng.ns ng.m]);
              for e = 1:ng.m
                edge = ng.Edges(e,:);
                if edge(1) < edge(2)
                  twin = [edge(2) edge(1)];
                  twinInd = ismember(ng.Edges,twin,'rows');
                  edgePO = rand([ng.ns ng.ns]);
                  ng.edgePayoffs(:,:,e) = edgePO;
                  ng.edgePayoffs(:,:,twinInd) = edgePO;
                end
              end
              ng.edgePayoffs = rand([ng.ns ng.ns ng.m]);
            case 'split'
              ng.edgePayoffs = zeros([ng.ns ng.ns ng.m]);
              numC = ceil(ng.m*ng.coordFraction);
              numAC = ng.m - numC;
              edgesC = randperm(ng.m,numC);
              edgesAC = setdiff(1:ng.m,edgesC);
              baseGameC = eye(ng.ns);
              baseGameC(1,1) = ng.coordLevel;
              baseGameAC = ones(ng.ns)-eye(ng.ns);
              baseGameAC(1,2) = ng.coordLevel;
              ng.edgePayoffs(:,:,edgesC) = repmat(baseGameC,[1 1 numC]) + ng.payoffVariance*rand([ng.ns ng.ns numC]);
              ng.edgePayoffs(:,:,edgesAC) = repmat(baseGameAC,[1 1 numAC]) + ng.payoffVariance*rand([ng.ns ng.ns numAC]);
              % epsEdges = rand(ng.m,1) < 0.5;
              % ng.edgePayoffs(:,:,epsEdges) = ng.edgePayoffs(:,:,epsEdges)*1e-6;
            case 'row-coord'
              ng.edgePayoffs = [2-cumsum(rand(1,2,ng.m));cumsum(rand(1,2,ng.m))];
          end
          ng.ns = size(ng.edgePayoffs(:,:,1),1);
        end
      end
      
      % Set update rules
      if isempty(ng.updateRules)
        switch ng.updateRule
          case 'Mixed'
            ng.updateRules = repmat({'BestResponse'},ng.n,1);
            randIm = rand(ng.n,1)>0.5;
            ng.updateRules(randIm) = {'Imitation'};
          otherwise 
            ng.updateRules = repmat({ng.updateRule},ng.n,1);
        end
      end
      
      % Set initial strategies
      switch ng.initType
        case 'uniform'
          ng.x0 = ng.ns*ones(ng.n,1);
        case 'random'
          rng(ng.seed+1);
          ng.x0 = ceil(ng.ns*rand(ng.n,1));
        case 'randomEq'
          rng(ng.seed+1);
          good = false;
          maxTries = 500;
          counter = 0;
          while ~good && counter < maxTries
            ng.x0 = ceil(ng.ns*rand(ng.n,1));
            xf = Play(ng);
            good = any(diff(xf)) && any(xf==1);
            counter = counter + 1;
          end
          ng.x0 = xf;
      end
      
      % Set payoff advantages (delta values) 2x2 games only
      deltaA = zeros(ng.m,1);
      deltaB = zeros(ng.m,1);
      
      if ng.ns==2
        for e = 1:size(ng.Edges,1)
          pi_ij = ng.edgePayoffs(:,:,e);
          deltaA(e) = pi_ij(1,1) - pi_ij(2,1);
          deltaB(e) = pi_ij(2,2) - pi_ij(1,2);
        end
      end
      
      % Store neighbor payoffs for faster updates
      [I,J] = find(ng.Adj);
      ng.deltaA = sparse(I,J,deltaA,ng.n,ng.n);
      ng.deltaB = sparse(I,J,deltaB,ng.n,ng.n);
      ng.Edges = [I J];
      ng.edgeAdjMap = sparse(I,J,1:length(I),ng.n,ng.n);
      ng.weightedAdj = sparse(ng.Adj.*dist(ng.Vxy'));
      
      ng.x = ng.x0;
      ng.deg = sum(ng.Adj,2);
      if ng.showGames
        InitPlots(ng,1);
      end
      ComputePayoffs(ng);
      % ComputeEnergy(ng);
    end
    
    function Draw(ng)
      if isempty(ng.hAgents)
        InitPlots(ng,1);
      end
      set(ng.hAgents,'CData',ng.ns-ng.x');
      if ng.nodeLabels
        for i = 1:length(ng.hText)
          set(ng.hText(i),'String',num2str(ng.y(i),'%0.0f'));
        end
      end
      drawnow;
      pause(.01);
      if ng.pauseGames
        pause;
      end
    end
    
    function disp(ng)
      disp(ng.Adj);
      disp(ng.payoffs);
    end
    
    function new = copy(ng)
      % Instantiate new object of the same class.
      new = feval(class(ng));
      
      % Copy all non-hidden properties.
      p = properties(ng);
      for i = 1:length(p)
        new.(p{i}) = ng.(p{i});
      end
      new.showGames = false;
    end
    
    function AddPayoff(ng,rowInd,agents,incentive)
      
      ng.payoffs(rowInd,:,agents) = ng.payoffs(rowInd,:,agents) + incentive;
      
      % Update edge payoffs
      ng.edgePayoffs = zeros(ng.ns,ng.ns,ng.m);
      for e = 1:size(ng.Edges,1)
        i = ng.Edges(e,1);
        pi_ij = ng.payoffs(:,:,i);
        ng.edgePayoffs(:,:,e) = pi_ij;
      end
    end
    
    function Reset(ng)
      ng.x = ng.x0;
      ComputePayoffs(ng);
    end
    
    function SaveFig(ng,str)
      if nargin < 2
        str = '';
      end
      set(gcf,'Color','w');
      export_fig([ng.saveStr str],'-pdf','-q101');
    end
    
    InitPlots(ng,fig,ax)
    [xf,terminated,phi] = Play(ng)
    ComputePayoffs(ng,agents)
    Phi = ComputeEnergy(ng,Adj,X,Xind,selNodes)
  end
end


