function A = hss_symmetrize(A, tol)
% Returns a version with symmetrized generators of the HSS matrix A 
if ~exist('tol', 'var')
	tol = hssoption('threshold')	;
end
if A.leafnode == 1
	return
end
	
if A.A11.leafnode == 1
	% A.A11.U = A.A11.V;
	% A.A22.V = A.A22.U;
	% A.B12 = A.B21';
	if norm(A.A11.D - A.A11.D', 'fro') > tol * norm(A.A11.D, 'fro') ||	norm(A.A22.D - A.A22.D', 'fro') > tol * norm(A.A22.D, 'fro')
		error('HSS_SYMMETRIZE:: trying to symmetrize a non Hermitian matrix\n');
	else
		A.A11.D = (A.A11.D + A.A11.D')/2;
		A.A22.D = (A.A22.D + A.A22.D')/2;
	end

	% X = A.A11.U \ A.A11.V;
	X = hss_aux_ls(A.A11.U, A.A11.V, tol);

	A.A11.U = A.A11.V;
	A.B12 = X * A.B12;

	if A.topnode ~= 1
		A.Rl = X * A.Rl;
	end

	% X = A.A22.V \ A.A22.U;
	X = hss_aux_ls(A.A22.V, A.A22.U, tol);
	A.A22.V = A.A22.U;
	A.B12 = A.B12 * X';

	if A.topnode ~= 1
		A.Wr = X * A.Wr;
	end

	if norm(A.B12 - A.B21', 'fro') > tol * norm(A.B21, 'fro')
		error('HSS_SYMMETRIZE:: trying to symmetrize a non Hermitian matrix\n');
	else
		A.B12 = A.B21';
	end	
else
	A.A11 = hss_symmetrize(A.A11, tol);
	A.A22 = hss_symmetrize(A.A22, tol);

	% X = [ A.A11.Rl ; A.A11.Rr ] \ [ A.A11.Wl ; A.A11.Wr ];
	X = hss_aux_ls([ A.A11.Rl ; A.A11.Rr ], [ A.A11.Wl ; A.A11.Wr ], tol);
	A.A11.Rl = A.A11.Wl; 
	A.A11.Rr = A.A11.Wr;
	A.B12 = X * A.B12;

	if A.topnode ~= 1
		A.Rl = X * A.Rl;
	end

	% X = [ A.A22.Wl ; A.A22.Wr ] \ [ A.A22.Rl ; A.A22.Rr ];
	X = hss_aux_ls([ A.A22.Wl ; A.A22.Wr ], [ A.A22.Rl ; A.A22.Rr ], tol);
	A.A22.Wl = A.A22.Rl; 
	A.A22.Wr = A.A22.Rr;
	A.B12 = A.B12 * X';

	if A.topnode ~= 1
		A.Wr = X * A.Wr;
	end

	if norm(A.B12 - A.B21', 'fro') > tol * norm(A.B21, 'fro')
		error('HSS_SYMMETRIZE:: trying to symmetrize a non Hermitian matrix\n');
	else
		A.B12 = A.B21';
	end	
end
end

function X = hss_aux_ls(A, B, tol)
	X = A \ B;

	if norm(A*X - B, 'fro') > tol * (norm(A, 'fro') + norm(B, 'fro'))
		error('HSS_SYMMETRIZE:: trying to symmetrize a non Hermitian matrix\n');
	end
end
