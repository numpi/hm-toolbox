function B = hss_build_diagonal(varargin)
%HSS_BUILD_DIAGONAL  Build a diagonal HSS matrix from a vector.
%
% B = HSS_BUILD_DIAGONAL(v) builds an HSS representation of diag(v).
%
switch nargin
    case 1
        v = varargin{1};
        if isvector(v)
            if length(v) < hssoption('block-size')
                B = hss(diag(v));
                return;
            end
            B = hss_build_diagonal_ric(v);
            B.topnode = 1;
        else
            error('HSS::Unsopported non vector argument')
        end
    otherwise
        error('Unsupported number of parameters');
end
end

function B = hss_build_diagonal_ric(v)
B = hss();
B.topnode = 0;
B.leafnode = 0;
m = length(v);
if m <= hssoption('block-size')
    B.leafnode = 1;
    B.D = diag(v);
    B.leafnode = 1;
    B.U = zeros(m, 0);
    B.V = zeros(m, 0);
else
    mmid = ceil(m/2);
    B.ml = mmid;
    B.nl = mmid;
    B.mr = m - mmid;
    B.nr = B.mr;
    B.A11 = hss_build_diagonal_ric(v(1:mmid));
    B.A22 = hss_build_diagonal_ric(v(mmid+1:end));
end

end
