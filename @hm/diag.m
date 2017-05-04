function d = diag(H)
%DIAG Obtain the diagonal of an H-matrix. 

if ~isempty(H.F)
	d = diag(H.F);
else
	d = [ diag(H.A11) ; diag(H.A22) ];
end


end

