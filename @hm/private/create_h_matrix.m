function obj = create_h_matrix(A)
%CREATE_H_MATRIX Given a dense matrix A, construct a hierarchical
%representation for it. 

min_block_size = hmoption('block-size');

obj = hm();

obj.F = [];
obj.sz = size(A);

if size(A, 1) <= min_block_size && size(A, 2) <= min_block_size
	obj.F = A;
else
	% Get the middle point
	mp = ceil(size(A, 1) / 2);
	
	obj.A11 = create_h_matrix(A(1:mp,1:mp));
	obj.A22 = create_h_matrix(A(mp+1:end,mp+1:end));
	
	[obj.U21, obj.V21] = compress_matrix(A(mp+1:end,1:mp));
	[obj.U12, obj.V12] = compress_matrix(A(1:mp,mp+1:end));
end


end

