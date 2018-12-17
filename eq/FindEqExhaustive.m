function eqGlobal = FindEqExhaustive(ng)
% Find equilbria by checking all possilble states

allStates = LogMat(ng.n,(1:ng.ns)');
numStates = size(allStates,1);

isEq = false(numStates,1);
for i = 1:numStates
  display(['Checking state ' num2str(i) ' of ' num2str(numStates)]);
  x0 = allStates(i,:);
  ng.x0 = x0';
  xf = Play(ng);
  if all(x0' == xf)
    isEq(i) = true;
  end  
end

eqGlobal = allStates(isEq,:);