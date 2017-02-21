function obj = create_diagonal_h_matrix(D)
%CREATE_DIAGONAL_H_MATRIX Create an H-matrix with the specified diagonal. 

obj = hm();
n = length(D);

if length(D) <= hmoption('block-size')
	obj.F = diag(D);
	obj.sz = [ n, n ];
else
	mp = ceil(length(D) / 2);
	obj.A11 = create_diagonal_h_matrix(D(1:mp));
	obj.A22 = create_diagonal_h_matrix(D(mp+1:end));
	obj.U12 = zeros(mp, 0);
	obj.V12 = zeros(n - mp, 0);
	obj.U21 = zeros(n - mp, 0);
	obj.V21 = zeros(mp, 0);
	
	obj.sz = [ n, n ];
end


end

