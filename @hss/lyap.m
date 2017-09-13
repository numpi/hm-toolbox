function X = lyap(A, C, varargin)
	if isempty(varargin)
		X = hss_dac_lyap(A, A, C);
	else
		X = hss_sparse_dac_lyap(A,A,C,varargin{1},varargin{1});
	end

end
