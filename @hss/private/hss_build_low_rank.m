function B = hss_build_low_rank(obj, varargin)
% Build the HSS representantion of a low rank matrix in the form U*V' or U*S*V'
U = varargin{1};
V = varargin{2};
B = obj;

if size(U,2) ~= size(V,2)
    error('HSS_BUILD_LOW_RANK:: dimensions of U and V not compatible')
end

if nargin == 4
    S = varargin{3};
    if size(S,1) ~= size(S,2) || size(S,1)~=size(U,2)
        error('HSS_BUILD_LOW_RANK:: dimensions of S not compatible')
    end
else
    S = eye(size(U,2));
end

if B.leafnode
    B.D = U * V';
    return;
end

B = hss_build_low_rank_ric(B, U, V, S);
B.topnode = 1;

end

function B = hss_build_low_rank_ric(B, U, V, S)

B.topnode = 0;
if B.leafnode
    B.U = U;
    B.V = V;
    B.D = U * S * V';
else
    mmid = B.ml;
    nmid = B.nl;
    
    B.B12 = S;
    B.B21 = S;
    B.Rl = eye(size(S,1)); B.Rr = eye(size(S,1));
    B.Wl = eye(size(S,1)); B.Wr = eye(size(S,1));
    
    B.A11 = hss_build_low_rank_ric(B.A11, U(1:mmid,:) , V(1:nmid,:), S);
    B.A22 = hss_build_low_rank_ric(B.A22, U(mmid+1:end,:), V(nmid+1:end,:), S);
end
end
