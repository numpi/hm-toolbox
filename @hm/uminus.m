function H = uminus(H1)
%UMINUS Change the sign of a HODLR Matrix 
if ~isempty(H1.F)
		H = H1;
		H.F = -H1.F;
else
	H = H1;
	H.A11 = -H1.A11;
	H.A22 = -H1.A22;
	H.U21 = -H1.U21;
	H.U12 = -H1.U12;
end
