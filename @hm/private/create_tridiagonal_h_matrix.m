function H = create_tridiagonal_h_matrix(A)
%CREATE_TRIDIAGONAL_H_MATRIX Create a tridiagonal H-matrix

H = hm();

block_size = hmoption('block-size');

if size(A, 1) <= block_size
	H.F = A;
	H.sz = size(A);
else
	mp = ceil(size(A, 1) / 2);
	n = size(A, 1);
	
	H.A11 = create_tridiagonal_h_matrix(A(1:mp,1:mp));
	H.A22 = create_tridiagonal_h_matrix(A(mp+1:end,mp+1:end));
	
	H.U12 = [ zeros(mp-1,1) ; A(mp,mp+1) ];
	H.V12 = [ 1 ; zeros(n - mp - 1, 1) ];
	
	H.U21 = [ A(mp+1,mp) ; zeros(n - mp - 1, 1) ];
	H.V21 = [ zeros(mp-1,1) ; 1 ];
	
	H.sz = size(A);
end


end

