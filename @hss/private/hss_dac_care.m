function X = hss_dac_care(A, B, C, spA, tol, debug, nrmtype)
% solve the CARE A'X + X A - X B B' X + C = 0 with a divide and conquer method
%
%  Assumptions:
%  A   hss matrix
%  spA sparse version of A (if any)
%  B   dense block column vector
%  C   Hermitian hss matrix (same clustering as A)
%------------------------------------------------------------------------------
if ~exist('debug', 'var') || isempty(debug)
	debug = false;
end
if ~exist('tol', 'var') || isempty(tol)
	tol = hssoption('threshold') * 1e2;
end
if ~exist('nrmtype', 'var') || isempty(nrmtype)
	nrmtype = 2;
end
if ~exist('spA','var') 
    spA = [];
end
[n, m] = size(A);
X = hss();
if  ~isempty(A.D) || ~isempty(C.D)
	X.leafnode = 1;
	X.topnode = 1;
	if ~exist('icare', 'file')
		X.D = care(A.D, B, C.D);
	else
		X.D = small_care_solve(A.D, B, C.D, hssoption('threshold'), 50);
		%X.D = icare(A.D, B, C.D);
	end

	if debug
        	fprintf('Base case: Dimension = %d, Residual = %e\n', size(A, 1), norm(A' * X + X * A - X * B * B' * X + C )/norm(X));
    	end
else
	A.topnode = 1;
	C.topnode = 1;
	nn = size(A.A11, 1);
	% Solve the block diagonal equation via recursion
	if issparse(spA)
        	X = blkdiag( hss_dac_care(A.A11, B(1:nn, :), C.A11, spA(1:nn, 1:nn),tol, debug, nrmtype),...
                     hss_dac_care(A.A22, B(nn + 1:end, :), C.A22, spA(nn + 1:end, nn + 1:end), tol, debug, nrmtype) );
	else
		X = blkdiag( hss_dac_care(A.A11, B(1:nn, :), C.A11, [],tol, debug, nrmtype),...
                     hss_dac_care(A.A22, B(nn + 1:end, :), C.A22, [], tol, debug, nrmtype) );
	end

        % Build low-rank representation of the correction equation
    	[CV, CU] = hss_offdiag(C, 'upper'); 
    	dC = kron([0, 1; 1, 0], eye(size(CU, 2)));
    	BU = blkdiag(B(1:size(A.A11, 1), :), B(size(A.A11, 1) + 1:end, :));
        [AU12, AV12] = hss_offdiag(A, 'upper');
        [AU21, AV21] = hss_offdiag(A, 'lower');
    	AU = [zeros(size(AU12, 1), size(AU21, 2)), AU12; AU21, zeros(size(AU21, 1), size(AU12, 2))];
    	AV = blkdiag(AV21, AV12);

        u = [blkdiag(CV, CU), AV, X * AU, X * BU];
        D = blkdiag(dC, kron([0, 1; 1, 0], eye(size(AU, 2))), -kron([0, 1; 1, 0], eye(size(BU, 2)/2)));
    	[u, D] = sym_compress_factors(u, D, hssoption('threshold'), 1); % re-compression step
    
    	if  ~isempty(u) && ~isempty(D)
		XB = -X' * B;
		if issparse(spA) % Build the struct for A' exploiting sparsity + low-rank structure
			[LA, UA, pA, qA] = lu(spA, 'vector');
        		iqA(qA) = 1:size(spA,1);
        		ipA(pA) = 1:size(spA,1);

			ATstruct = struct(...
        	        'solve', @(nu, mu, x) sparse_woodbury_solve(nu, mu, UA', LA', qA, ipA, x, XB, B), ...
        	        'multiply', @(rho, eta, x) rho * spA' * x + rho * XB * (B' * x) - eta * x, ...
        	        'isreal', isreal(spA), ...
        	        'nrm', normest(spA, 1e-2));
		else
			A.topnode = 1;
			ATstruct = ek_struct(A' + hss('low-rank', XB, B)); % Build the structure exploiting the HSS structure
		end
		[dXU, ~, dUZ,~, rkres] = ek_care(ATstruct, B, u, min(1000, n), tol, debug, nrmtype, 'kernel', D); % Solve the CARE with low-rank rhs
    	else
        	dXU = []; dUZ = [];
    	end
    	if ~isempty(dXU) && ~isempty(dUZ)
        	if debug
            		dX = hss('low-rank', dXU * dUZ, dXU); F = hss('low-rank', B, B);
            		Ac = A - F * X; Q = u * D * u';
            		fprintf('Correction eq: Dimension = %d, Rank of the sol. = %d, Residual = %e\n', size(A, 1), size(dXU, 2), norm(Ac' * dX + dX * Ac - dX * F * dX + Q) / norm(dX))
        	end
		X = X + hss('low-rank', dXU * dUZ, dXU);
		X.topnode = A.topnode;
		if debug 
        		fprintf('Complete equation: Dimension = %d, HSS rank sol. = %d, Residual = %e\n', size(A, 1), hssrank(X), norm(A' * X + X * A - X * F * X + C, 1e-2) / norm(X));
    		end
    	end
end
end

function [u, D] = sym_compress_factors(u, D, tol, nrm)
    [qu, ru] = qr(u, 0);
    [uu, s] = eig(ru * D * ru');
    [~, ord] = sort(diag(abs(s)), 'descend');
    s = s(ord, ord); uu = uu(:, ord);
    rk = sum(diag(abs(s)) > tol * abs(s(1,1) * nrm));
    u = qu * uu(:, 1:rk);
    D = s(1:rk, 1:rk);
end


