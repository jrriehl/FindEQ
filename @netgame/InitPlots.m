function InitPlots(netgame,fig,ax)
% Initialize plot of network game

n = netgame.n;
ns = netgame.ns;
x0 = netgame.x0;
normInd = (1:n)';
E = netgame.Edges;
Vxy = netgame.Vxy;
dim = size(Vxy,2);

if ~isfield(netgame,'nodeSize')
  netgame.nodeSize = 10;
end

scatterSize = pi*netgame.nodeSize^2/4;
markerSize = netgame.nodeSize;
markerEdgeWidth = markerSize/8;

if nargin == 2
  figure(fig); clf; hold on;
elseif nargin == 1
  figure(1); clf; hold on;
else
  axes(ax); hold on;
end

meanEdgeDir = zeros(n,dim);
for e = 1:size(E,1)
  i = E(e,1); j = E(e,2);
  vi = Vxy(i,:);
  vj = Vxy(j,:);
  if dim == 2
    line([vi(1) vj(1)],[vi(2) vj(2)],'Color',[.5 .5 .5],'LineWidth',1.25);
    % arrow([vi(1) vi(2)],[vj(1) vj(2)],'Color',[.5 .5 .5]);
  elseif dim == 3
    line([vi(1) vj(1)],[vi(2) vj(2)],[vi(3) vj(3)],'Color',[.5 .5 .5],'LineWidth',1.25);
    % arrow([vi(1) vi(2)],[vj(1) vj(2)],'Color',[.5 .5 .5]);
  end
  meanEdgeDir(i,:) = meanEdgeDir(i,:) + (vj-vi)/norm(vj-vi);
  meanEdgeDir(j,:) = meanEdgeDir(j,:) + (vi-vj)/norm(vi-vj);
end
meanEdgeDir = meanEdgeDir./repmat(sqrt(sum(meanEdgeDir.^2,2)),1,dim);
meanEdgeDir(isnan(meanEdgeDir)) = 0;

if dim == 2
  netgame.hAgents = scatter(Vxy(:,1),Vxy(:,2),scatterSize,ns-x0,'filled','Marker','o');
  netgame.hNormal = plot(Vxy(normInd,1),Vxy(normInd,2),'o','Color',[.4 .4 .4],'MarkerSize',markerSize,'LineWidth',markerEdgeWidth);
  
  % Init. active and switching agents
  netgame.hActive = plot(0,0,'bo','MarkerSize',netgame.nodeSize,'LineWidth',2.5);
  set(netgame.hActive,'XData',[],'YData',[]);
  netgame.hSwitch = plot(0,0,'o','Color',[0 .7 0],'MarkerSize',netgame.nodeSize,'LineWidth',2.5);
  set(netgame.hSwitch,'XData',[],'YData',[]);
  
  % Initialize agent text
  if netgame.nodeLabels
    dx = .04*meanEdgeDir(:,1)+.02;
    dy = .04*meanEdgeDir(:,2);
    netgame.hText = text(Vxy(:,1)-dx,Vxy(:,2)-dy,num2str((1:n)'),'FontSize',16);
    % netgame.hText = text(Vxy(:,1)-dx,Vxy(:,2)-dy,num2str((1:n)'),'FontSize',16);
  end
  
elseif dim == 3
  netgame.hAgents = scatter3(Vxy(:,1),Vxy(:,2),Vxy(:,3),scatterSize,ns-x0,'filled','Marker','o');
  netgame.hNormal = plot3(Vxy(normInd,1),Vxy(normInd,2),Vxy(normInd,3),'o','Color',[.4 .4 .4],'MarkerSize',markerSize,'LineWidth',markerEdgeWidth);

  if netgame.nodeLabels
    dx = .04*meanEdgeDir(:,1) + .02;
    dy = .04*meanEdgeDir(:,2) - .001;
    dz = .04*meanEdgeDir(:,3) - .001;
    netgame.hText = text(Vxy(:,1)-dx,Vxy(:,2)-dy,Vxy(:,3)-dz,num2str((1:n)'),'FontSize',16);
    % netgame.hText = text(Vxy(:,1)-dx,Vxy(:,2)-dy,num2str((1:n)'),'FontSize',16);
  end

end

caxis([0 ns-1]);
colormap(gray);
axis equal;
axis off;

if netgame.saveFigures
  SaveFig(netgame,'-Init')
end

end