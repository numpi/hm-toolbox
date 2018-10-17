function l = eigs(A, eig_type, tol)
	if ~exist('tol', 'var')
		tol = 1e-6;
	end
	if ~exist('eig_type', 'var')
		eig_type = 'LM';
	end

	if strcmp(eig_type, 'SM') == 1
		l = hss_inv_pow_met(A, tol);
	elseif strcmp(eig_type, 'LM') == 1
		l = hss_pow_met(A, tol);
	else
		error('EIGS:: not supported eig_type')
	end
	

