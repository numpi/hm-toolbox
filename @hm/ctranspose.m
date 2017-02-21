function Ht = ctranspose(H)
%TRANSPOSE Conjugate transpose of the H-matrix H. 

Ht = H;

if ~isempty(H.F)
	Ht.F = H.F';
else
	Ht.A11 = H.A11';
	Ht.A22 = H.A22';
	Ht.U12 = H.V21;
	Ht.V12 = H.U21;
	Ht.U21 = H.V12;
	Ht.V21 = H.U12;
end

end
