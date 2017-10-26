function X = lyap(A, C, varargin)
	if isempty(varargin)
		X = hss_dac_lyap(A, A, C);
	else
		nrmA = 1 / norm(A, 'fro');
		
		A = A * nrmA;
		C = C * nrmA;
		sA = varargin{1} * nrmA;
		sB = varargin{2} * nrmA;
								
		X = hss_sparse_dac_lyap(A, A, C, sA, sB);
	end

end
