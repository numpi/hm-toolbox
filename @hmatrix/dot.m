function y = dot(A, B)
% Compute trace(A' * B) with either A or B HMATRIX

if ~isa('A', 'hmatrix')
	HA = hmatrix; HA.sz = size(A); HA.F = A;
	A = HA;
end	
if ~isa('B', 'hmatrix')
	HB = hmatrix; HB.sz = size(B); HB.F = B;
	B = HB;
end	

if A.admissible  % A is a low-rank leaf node
	V = B * A.V;
	y = trace(U' * V);
elseif B.admissible % B is a low-rank leaf node
	U = A' * B.U;
	y = trace(B.V' * U);
elseif is_leafnode(A) && is_leafnode(B) % Both are dense leaf nodes
	y = trace(A.F' * B.F);
elseif ~is_leafnode(A) && ~is_leafnode(B) % Both are HMATRIX
	y = dot(A.A11, B.A11) + dot(A.A21, B.A21) + dot(A.A12, B.A12) + dot(A.A22 * B.A22);
else % One of the two is a dense leaf node and the other is an HMATRIX
	if is_leafnode(A)
		[m1, n1] = size(B.A11);
		y = dot(A.F(1:m1, 1:n1), B.A11) + dot(A.F(m1+1:end, 1:n1), B.A21) + dot(A.F(1:m1, n1+1:end), B.A12) + dot(A.F(m1+1:end, n1+1:end), B.A22);
	else
		[m1, n1] = size(A.A11);
		y = dot(A.A11, B.F(1:m1, 1:n1)) + dot(A.A21, B.F(m1+1:end, 1:n1)) + dot(A.A12, B.F(1:m1, n1+1:end)) + dot(A.A22, B.F(m1+1:end, n1+1:end));
	end
end
