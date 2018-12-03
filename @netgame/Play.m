function [x,t,phi] = Play(ng)
% Simulate the network game
ng.x = ng.x0;
ComputePayoffs(ng);

baseSaveStr = ng.saveStr;

% Draw initial strategies
if ng.showGames
  Draw(ng);
  if ng.saveFigures
    SaveFig(ng,'-00');
  end
end

phi = [];
t = 0;
ng.terminated = false;
while ~ng.terminated && t < ng.numRounds
  
  ng.saveStr = [baseSaveStr '-' num2str(t,'%0.2d')];
  
  % Update strategies
  Update(ng);
  
  % Compute current round payoffs
  ComputePayoffs(ng);
  
  t = t + 1;
  
  if ng.showEnergy
    phit = ComputeEnergy(ng);
    phi = [phi;phit];
    display(['Energy: ' num2str(sum(ng.phi))]);
  end
  
  if ng.showGames
    Draw(ng);
    if ng.saveFigures
      SaveFig(ng,'-PostSwitch');
    end
    if ng.showActivations
      set(ng.hActive,'XData',[],'YData',[]);
      set(ng.hSwitch,'XData',[],'YData',[]);
    end
  end
end

if ~ng.terminated
  disp(['Game did not terminate after ' num2str(t) ' rounds.']);
  t = inf;
else
  if ng.saveFigures
    SaveFig(ng,'-Final');
  end
end

x = ng.x;
ng.saveStr = baseSaveStr;