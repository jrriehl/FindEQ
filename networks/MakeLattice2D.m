function [A,Vxy] = MakeLattice2D(nrows,ncols)

n = nrows*ncols;
V = (1:n)';
A = zeros(n,n);

Vright = V(rem(V,ncols)==0);
Vleft = V(rem(V,ncols)==1);
Vtop = V(1:ncols);
Vbottom = V(n-ncols+1:n);
Vcorners = unique([1;ncols;n-ncols+1;n]);
Vleft = setdiff(Vleft,Vcorners);
Vright = setdiff(Vright,Vcorners);
Vtop = setdiff(Vtop,Vcorners);
Vbottom = setdiff(Vbottom,Vcorners);
Vint = setdiff(V,[Vright;Vleft;Vtop;Vbottom;Vcorners]);

Eint = [[Vint Vint+1];[Vint-1 Vint];[Vint Vint+ncols];[Vint Vint-ncols]];
Eleft = [[Vleft Vleft+1];[Vleft Vleft+ncols];[Vleft Vleft-ncols]];
Eright = [[Vright-1 Vright];[Vright Vright+ncols];[Vright Vright-ncols]];
Etop = [[Vtop Vtop+1];[Vtop-1 Vtop]];
Ebottom = [[Vbottom Vbottom+1];[Vbottom-1 Vbottom]];
Ecorners = [[1 2];[n-1 n]];
if nrows > 1 && ncols > 1
  Etop = [Etop;[Vtop Vtop+ncols]];
  Ebottom = [Ebottom;[Vbottom Vbottom-ncols]];
  Ecorners = [Ecorners; [1 ncols+1];[ncols-1 ncols];[ncols 2*ncols];...
    [n-2*ncols+1 n-ncols+1];[n-ncols+1 n-ncols+2];[n-ncols n]];
end  

E = unique([Eint;Eleft;Eright;Etop;Ebottom;Ecorners],'rows');

for e = 1:size(E,1)
  i = E(e,1);
  j = E(e,2);
  A(i,j) = 1;
  A(j,i) = 1;
end

qr = (1/nrows:1/nrows:1) - 1/(2*nrows);
qc = (1/ncols:1/ncols:1) - 1/(2*ncols);
Vx = kron(qr,ones(ncols,1));
Vy = kron(qc,ones(nrows,1))';
Vxy = [Vx(:) Vy(:)];

