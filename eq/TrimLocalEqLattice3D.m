function TrimLocalEqLattice3D(ng)

deg = ng.deg;
buried = deg == 6;

% Trim local equilibria
for i = 1:length(ng.LEQ)
  numLoc = size(ng.LEQ{i},1);
  keep = true(numLoc,1);
  for e = 1:size(ng.LEQ{i},1)
    % According to current energy function, buried nodes will always be 3
    eqInds = sort([i find(ng.Adj(i,:))]);
    xe = ng.LEQ{i}(e,:);
    if any(ismember(xe(buried(eqInds)),[1 2]))
      keep(e) = false;
    end
    if any(xe(~buried(eqInds)) == 3)
      keep(e) = false;
    end
  end
  ng.LEQ{i} = ng.LEQ{i}(keep,:);
end