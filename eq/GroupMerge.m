function [MergedREQ,Vgroup] = GroupMerge(REQ,selInds,PiPart)
% Merge a set of regional equilibrium sets such that the output set
% consists of all compatible combinations of the input sets

if nargin < 2
  selInds = 1:length(REQ);
end
numGroups = length(selInds);
n = size(PiPart,1);

MergedREQ = uint8(zeros(0,0));

Vgroup = false(n,0);

% Loop through all groups to check for compatibility of regionally stable
% strategy states
for ii = 1:numGroups
  i = selInds(ii);
  if isempty(MergedREQ)
    Vgroup = PiPart(:,ii);
    MergedREQ = REQ{i};
  elseif isempty(REQ{i})
    MergedREQ = uint8(zeros(0,0));
    return;
  else
    Vi = PiPart(:,ii); 
    MergedREQ = Merge(MergedREQ,REQ{i},Vgroup,Vi);
    Vgroup = Vgroup | Vi;
  end
end