function B = hss_from_banded(obj, varargin)
%HSS_FROM_BANDED Build an HSS matrix from a banded one.
%
% B = HSS_FROM_BANDED(A) builds an HSS representation of the banded matrix A.
%     The bandwidth is determined automatically. If possible, consider passing
%     arguments indicating the lower and upper bandwidth, in the form
%
%        B = HSS_FROM_BANDED(A, BL, BU);
%

switch nargin
    case 2
        A = varargin{1};
        [bl, bu] = bandwidth(A);
        
    case 4
        A = varargin{1};
        bl = varargin{2};
        bu = varargin{3};
        
    otherwise
        error('Unsupported number of parameters');
end

B = obj;

if B.leafnode
    B.D = full(A);
    return;
end

B = hss_build_band_ric(B, A, bl, bu);
B.topnode = 1;

end

function B = hss_build_band_ric(B, A, bl, bu)

if B.leafnode
    B.D = full(A);
    
    [ml, nl] = size(A);
    
    if ml ~= nl
        error(['The banded constructor is only supported for ' ...
            'matrices with square diagonal blocks; please use hss(A, ' ...
            char(39) 'cluster' char(39) ', rc, cc) instead' ]);
    end
    
    B.U = zeros(ml,bu+bl);
    B.U(1:bl,1:bl) = eye(bl);
    B.U(end-bu+1:end,bl+1:end) = eye(bu);
    
    B.V = zeros(nl,bu+bl);
    B.V(1:bu,1:bu) = eye(bu);
    B.V(end-bl+1:end,bu+1:end) = eye(bl);
else
    mmid = B.ml;
    nmid = B.nl;
    
    B.Rl = zeros(bl+bu); B.Rl(1:bl,1:bl) = eye(bl);
    B.Rr = zeros(bl+bu); B.Rr(end-bu+1:end,end-bu+1:end) = eye(bu);
    B.Wl = zeros(bl+bu); B.Wl(1:bu,1:bu) = eye(bu);
    B.Wr = zeros(bl+bu); B.Wr(end-bl+1:end,end-bl+1:end) = eye(bl);
    B.B21 = [zeros(bl,bu),full(A(mmid + 1:mmid + bl, nmid - bl + 1:nmid));zeros(bu,bl+bu)];
    B.B12 = [ zeros(bl, bl + bu);full(A(mmid-bu+1:mmid,nmid+1:nmid+bu)),zeros(bu,bl); ];
    
    B.A11 = hss_build_band_ric(B.A11, A(1:mmid,1:nmid), bl, bu);
    B.A22 = hss_build_band_ric(B.A22, A(mmid+1:end,nmid+1:end), bl, bu);
end
end
