function A = hm2hss(B)
%HM2HSS Conversion between data-sparse formats. 

m = size(B, 1);
n = size(B, 2);

A = hss('banded', sparse(m, n), 0, 0);

if hmoption('block-size') >= min(m, n)
	A = hss(full(B));
else
	UU = [ B.U12 , zeros(size(B.U12,1), size(B.U21,2)) ; ...
			  zeros(size(B.U21, 1), size(B.U12, 2)) , B.U21 ];
	VV = [ zeros(size(B.V12, 1), size(B.V21, 2)) , B.V21 ; ... 
			  B.V12 , zeros(size(B.V12,1), size(B.V12,2))];
	
	A = blkdiag(hm2hss(B.A11), hm2hss(B.A22)) + ...
			hss('low-rank', UU, VV);
end

%A = compress(A);

end

