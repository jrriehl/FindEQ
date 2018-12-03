function [hs,ht] = NetworkPlot(Vxy,Edges,values,ax)
% Vxy = 2-column array of X and Y coordinates for nodes
% Edges = 2-column array of edges [I J]
% values = node values to be displayed

n = size(Vxy,1);

if nargin < 4
  figure(1); clf; hold on;
  ax = gca;
end

for i = 1:length(Edges)
  vi = Vxy(Edges(i,1),:);
  vj = Vxy(Edges(i,2),:);
  line([vi(1) vj(1)],[vi(2) vj(2)],'Color',[.5 .5 .5]);
  hold on;
end
hs = scatter(ax,Vxy(:,1),Vxy(:,2),50,1-values,'filled','Marker','o');
plot(ax,Vxy(:,1),Vxy(:,2),'o','Color',[.4 .4 .4],'MarkerSize',8,'LineWidth',2);
% ht = text(Vxy(:,1)-.1,Vxy(:,2),num2str(defStrat,'%0.2f'));
colormap(gray);
axis equal;
axis off;

end

