% Build the HSS representantion of a low rank matrix in the form U*V' or U*S*V'
function B = hss_build_low_rank(varargin)
		U = varargin{1};
		V = varargin{2};
		B = hss();
		if size(U,2) ~= size(V,2)
			error('HSS_BUILD_LOW_RANK:: dimensions of U and V not compatible')
		end
		if nargin == 3
			S = varargin{3};
			if size(S,1) ~= size(S,2) || size(S,1)~=size(U,2)
				error('HSS_BUILD_LOW_RANK:: dimensions of S not compatible')
			end
		else
			S = eye(size(U,2));
		end

                if size(U,1) < hssoption('block-size')
                       B = hss(U*V');
                       return;
                end

		k = 1+ceil(log2(size(U,1)/hssoption('block-size')));	
		B = hss_build_low_rank_ric(U, V, S, k);
		B.topnode = 1;	
		%B = rmfield(B, {'Rl', 'Rr', 'Wr', 'Wl'});
end

function B = hss_build_low_rank_ric(U, V, S, k)
	B = hss();
	B.topnode = 0;
	if k == 1
		B.leafnode = 1;
		B.U = U;
		B.V = V;
		B.D = U * S * V';
	else
		B.leafnode = 0;
		mmid = ceil(size(U,1)/2);
		nmid = ceil(size(V,1)/2);
		B.ml = mmid; B.nl = nmid;
		B.mr = size(U,1) - mmid; B.nr = size(V,1) - nmid;
		B.Bu = S;
		B.Bl = S;
		B.Rl = eye(size(S,1)); B.Rr = eye(size(S,1));
		B.Wl = eye(size(S,1)); B.Wr = eye(size(S,1));
		B.hssl = hss_build_low_rank_ric( U(1:mmid,:) , V(1:nmid,:), S, k-1);
		B.hssr = hss_build_low_rank_ric(U(mmid+1:end,:), V(nmid+1:end,:), S, k-1);
	end
end
