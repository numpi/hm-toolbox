function rk = hmrank(H)
%QSRANK Obtain the maximum rank in the off-diagonal blocks of H. 
%
% RK = HSRANK(H) returns the maximum rank in the representation of the
% off-diagonal blocks. 

if ~isempty(H.F)
	rk = 0;
else
	rk = max([ hmrank(H.A11), hmrank(H.A22), ...
			   size(H.U12, 2), size(H.U21, 2) ]);
end

end

