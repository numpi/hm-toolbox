function B = hss_build_diagonal(obj, varargin)
%HSS_BUILD_DIAGONAL  Build a diagonal HSS matrix from a vector.
%
% B = HSS_BUILD_DIAGONAL(v) builds an HSS representation of diag(v).
%
switch nargin
    case 2
        v = varargin{1};
        B = obj;
        if isvector(v)
            if B.leafnode
                B.D = diag(v);
                return;
            end
            B = hss_build_diagonal_ric(B, v);
            B.topnode = 1;
        else
            error('HSS::Unsopported non vector argument')
        end
    otherwise
        error('Unsupported number of parameters');
end
end

function B = hss_build_diagonal_ric(B, v)

m = length(v);

if B.leafnode == 1
    B.D = diag(v);
    B.U = zeros(m, 0);
    B.V = zeros(m, 0);
else
    if B.ml ~= B.nl
        error('Diagonal constructor is only supported for square matrices');
    end
    
    mmid = B.ml;
    
    B.A11 = hss_build_diagonal_ric(B.A11, v(1:mmid));
    B.A22 = hss_build_diagonal_ric(B.A22, v(mmid+1:end));
end

end
