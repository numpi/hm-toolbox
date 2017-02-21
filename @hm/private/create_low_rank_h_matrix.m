function obj = create_low_rank_h_matrix(U, V)
%CREATE_LOW_RANK_H_MATRIX Create a low rank H matrix. 

global hm_block_size

if isempty(hm_block_size)
	hmoption('block-size');
end

obj = hm();

obj.sz = [ size(U, 1), size(V, 1) ];

if obj.sz(1) <= hm_block_size
	obj.F = U * V';
else
	mp = ceil(obj.sz(1) / 2);
	obj.A11 = create_low_rank_h_matrix(U(1:mp,:), V(1:mp,:));
	obj.A22 = create_low_rank_h_matrix(U(mp+1:end,:), V(mp+1:end,:));
	obj.U12 = U(1:mp,:);
	obj.V12 = V(mp+1:end,:);
	obj.U21 = U(mp+1:end,:);
	obj.V21 = V(1:mp,:);
end


end

