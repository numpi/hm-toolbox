function rk = qsrank(H, exact)
%QSRANK Obtain the quasiseparability rank of H. 
%
% RK = QSRANK(H) returns the maximum rank in the representation of the
% off-diagonal blocks. 

if ~isempty(H.F)
	rk = 0;
else
	rk = max([ qsrank(H.A11), qsrank(H.A22), ...
			   size(H.U12, 2), size(H.U21, 2) ]);
end

end

