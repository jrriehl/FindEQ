function Vxy = trimtreeplot(p,c,d)
% TRIMTREEPLOT A modified standard function TREEPLOT. Plots the
%   leaves differently from TREEPLOT. They appear in their
%   respective layer instead of the deepest layer so that the tree
%   appears as "trimmed". The format is the same as for TREEPLOT,
%   see below. 
%
% trimtreeplot(p [,c,d])
%
% TREEPLOT Plot picture of tree.
%   TREEPLOT(p) plots a picture of a tree given a vector of
%   parent pointers, with p(i) == 0 for a root. 
%
%   TREEPLOT(P,nodeSpec,edgeSpec) allows optional parameters nodeSpec
%   and edgeSpec to set the node or edge color, marker, and linestyle.
%   Use '' to omit one or both.
%
%   See also ETREE, TREELAYOUT, ETREEPLOT.

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 5.8 $  $Date: 1997/11/21 23:44:59 $
%
%   Modified by RG 01-Nov-1 to use trimtreelayout instead of
%   TREELAYOUT 

% RG use another treelayout 
%[x,y,h]=treelayout(p);
[x,y,h]=trimtreelayout(p);

f = find(p~=0);
pp = p(f);
X = [x(f); x(pp); repmat(NaN,size(f))];
Y = [y(f); y(pp); repmat(NaN,size(f))];
X = X(:);
Y = Y(:);

if nargin == 1,
    n = length(p);
    if n < 500,
        plot (x, y, 'ro', X, Y, 'r-');
    else,
        plot (X, Y, 'r-');
    end;
else,
    [ignore, clen] = size(c);
    if nargin < 3, 
        if clen > 1, 
            d = [c(1:clen-1) '-']; 
        else,
            d = 'r-';
        end;
    end;
    [ignore, dlen] = size(d);
    if clen>0 & dlen>0
        plot (x, y, c, X, Y, d);
    elseif clen>0,
        plot (x, y, c);
    elseif dlen>0,
        plot (X, Y, d);
    else
    end;
end;
xlabel(['height = ' int2str(h)]);
axis([0 1 0 1]);

Vxy = [x' y'];