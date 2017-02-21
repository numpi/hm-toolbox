function H = create_banded_h_matrix(A, band)
%CREATE_BANDED_H_MATRIX Create a banded H-matrix. 

H = hm();

block_size = hmoption('block-size');

if size(A, 1) <= block_size
	H.F = A;
	H.sz = size(A);
else
	mp = ceil(size(A, 1) / 2);
	n = size(A, 1);
	
	if band <= min(n - mp)
		H.A11 = create_banded_h_matrix(A(1:mp,1:mp), band);
		H.A22 = create_banded_h_matrix(A(mp+1:end,mp+1:end), band);

		H.U12 = [ zeros(mp-band,band) ; A(mp-band+1:mp,mp+1:mp+band) ];
		H.V12 = [ eye(band) ; zeros(n - mp - band, band) ];

		H.U21 = [ A(mp+1:mp+band,mp-band+1:mp) ; zeros(n - mp - band, band) ];
		H.V21 = [ zeros(mp-band, band) ; eye(band) ];
	else
		H = create_h_matrix(A);
	end
	
	H.sz = size(A);
end


end

