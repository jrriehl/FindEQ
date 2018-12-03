function Update(ng,agents)

if nargin < 2
  agents = 1:ng.n;
end

xnew = ng.x;
xpm = 2*(2 - ng.x) - 1;

for ii = 1:length(agents)
  i = agents(ii);
  xi = ng.x(i);
  neighbors = find(ng.Adj(i,:))';
  switch ng.updateRules{i}
    
    case 'Imitation'
      [yMax,~] = max(ng.y(neighbors));
      if yMax > ng.y(i) % && ~ismember(i,ng.controlNodes),
        bestNeighbors = neighbors(ng.y(neighbors) == yMax);
        % If my strategy is not used by one of my best neighbors, I'll switch
        if ~any(bsxfun(@eq,xi,ng.x(bestNeighbors)))
          % randChoice = bestNeighbors(ceil(length(bestNeighbors)*rand));
          xnew(i) = ng.x(bestNeighbors(1));
        end
      end
      
    case 'FastIM'
      
      selfNeighbors = [i;neighbors];
      [~,j] = max(ng.y(selfNeighbors));
      xnew(i) = ng.x(selfNeighbors(j));
      
    case 'FastBR'
      
      hTest = ng.deltaA(i,neighbors)*(1+xpm(neighbors)) - ng.deltaB(i,neighbors)*(1-xpm(neighbors));
      if hTest ~= 0
        if hTest > 0
          xnew(i) = 1;
        else
          xnew(i) = 2;
        end
      end
      
    case 'BestResponse'
      % Compute payoffs of all strategies against current neighbor strategies
      
      % Find out-edges
      outEdgeInds = ng.edgeAdjMap(i,:);
      outNeighbors = find(outEdgeInds);
      outEdgeInds = outEdgeInds(outNeighbors);
      
      yTest = zeros(ng.ns,1);
      for k = 1:ng.ns        
        for j = 1:length(outNeighbors)
          out = outNeighbors(j);
          e = outEdgeInds(j);
          yTest(k) = yTest(k) + ng.edgePayoffs(k,ng.x(out),e);
        end
      end
      
      % Best response strategy
      [~,xbr] = max(yTest);
      
      % Only switch strategies if there is a strictly better response strategy
      if yTest(xbr) > yTest(ng.x(i))
        xnew(i) = xbr;
      end
      
    case 'BestResponseControl'
      
      % Find out-edges
      outEdgeInds = ng.edgeAdjMap(i,:);
      outNeighbors = find(outEdgeInds);
      outEdgeInds = outEdgeInds(outNeighbors);
      
      yTest = zeros(ng.ns,1);
      for k = 1:ng.ns
        for j = 1:length(outNeighbors)
          out = outNeighbors(j);
          e = outEdgeInds(j);
          yTest(k) = yTest(k) + ng.edgePayoffs(k,ng.x(out),e);
        end
      end
      
      % Best response strategy
      [~,xbr] = max(yTest);
      
      % Only switch strategies if there is a strictly better response
      % strategy and this node is not controlled
      if yTest(xbr) > yTest(ng.x(i)) && ~(ng.x(i)==3)
        xnew(i) = xbr;
      end
      
    case 'LogitResponse'
      % Soft threshold-based update
      
      % Find in-edges
      outEdges = ng.Edges(ng.Edges(:,1)==i,:);
      
      yTest = zeros(ng.ns,1);
      for k = 1:ng.ns
        for j = 1:size(outEdges,1)
          out = outEdges(j,1);
          yTest(k) = yTest(k) + ng.edgePayoffs(k,ng.x(out),i);
        end
      end
      
      % Compute logit choice utilities
      expUtils = exp(ng.logitBeta*yTest);
      normExpUtils = expUtils/sum(expUtils);
      cumSumNormExpUtils = cumsum(normExpUtils);
      
      % Select weighted random strategy
      r = rand;
      test = r < cumSumNormExpUtils;
      xnew(i) = find(test,1,'first');
      
    case 'ProteinBR'
      
      % Find out-edges
      outEdgeInds = ng.edgeAdjMap(i,:);
      outNeighbors = find(outEdgeInds);
      neighborStates = ng.x(outNeighbors);
      max_neighbors = 2*size(ng.Vxy,2);
      
      yTest = zeros(ng.ns,1);
      for k = 1:ng.ns
        yTest(k) = LocalProteinEnergy(k,neighborStates,max_neighbors);
      end
      
      %  Best response strategy
      [~,xbr] = min(yTest);
      
      % Only switch strategies if there is a strictly better response strategy
      if yTest(xbr) < yTest(ng.x(i))
        xnew(i) = xbr;
      end
      
    case 'Bank'
      
      % Find in-edges
      inEdgeInds = ng.edgeAdjMap(:,i);
      inNeighbors = find(inEdgeInds);
      neighborStates = ng.x(inNeighbors);
      
      pBar = ng.lambda;         % amount of full payments
      bankAlpha = ng.deltaA;    % pBar - fixed income
      bankBeta = ng.deltaB;     %
      
      incomes = (2-neighborStates) .* pBar(inNeighbors);
      
      if bankAlpha(i) <= bankBeta(inNeighbors,i)' * incomes
        xnew(i) = 1;
      else
        xnew(i) = 2;
      end
  end
end

active = xnew ~= ng.x;

if ~any(active)
  ng.terminated = true;  
else
  numActive = sum(active);
  activeInds = find(active);
  
  if ng.synchronous
    switchAgents = activeInds;
  else
    switchAgents = activeInds(ceil(rand*numActive));
  end
  
  if ng.showGames
    if ng.showActivations
      set(ng.hActive,'XData',ng.Vxy(active,1),'YData',ng.Vxy(active,2));
      if ng.pauseGames
        pause;
      end
      if ng.saveFigures
        SaveFig(ng,'-Active');
      end
      set(ng.hSwitch,'XData',ng.Vxy(switchAgents,1),'YData',ng.Vxy(switchAgents,2));
      if ng.pauseGames
        pause;
      end
      if ng.saveFigures
        SaveFig(ng,'-Switch');
      end
    end
  end
    
  ng.x(switchAgents) = xnew(switchAgents);
  
end

end
