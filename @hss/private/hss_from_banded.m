function B = hss_from_banded(varargin)
%HSS_FROM_BANDED Build an HSS matrix from a banded one.
%
% B = HSS_FROM_BANDED(A) builds an HSS representation of the banded matrix A.
%     The bandwidth is determined automatically. If possible, consider passing
%     arguments indicating the lower and upper bandwidth, in the form
%
%        B = HSS_FROM_BANDED(A, BL, BU);
%

switch nargin
    case 1
        A = varargin{1};
        %A
        [bl, bu] = bandwidth(A);
        
    case 3
        A = varargin{1};
        bl = varargin{2};
        bu = varargin{3};
        
    case 4
        A = varargin{2};
        bl = varargin{3};
        bu = varargin{4};
        
    otherwise
        error('Unsupported number of parameters');
end



if size(A,1) < hssoption('block-size')
    B = hss(A);
    return;
end

k = max(0, 1 + ceil(log2(size(A,1)/hssoption('block-size'))));
if k > 0
    B = hss_build_band_ric(A,bl,bu,k);
else
    B = hss();
    B.D = A;
end
B.topnode = 1;
% B = rmfield(B, {'Rl', 'Rr', 'Wr', 'Wl'});
end
function B = hss_build_band_ric(A,bl,bu,k)
B = hss();
B.topnode = 0;
B.leafnode = 0;
if k == 1
    B.leafnode = 1;
    B.D = full(A);
    
    [ml nl] = size(A);
    
    B.U = zeros(ml,bu+bl);
    B.U(1:bl,1:bl) = eye(bl);
    B.U(end-bu+1:end,bl+1:end) = eye(bu);
    
    B.V = zeros(nl,bu+bl);
    B.V(1:bu,1:bu) = eye(bu);
    B.V(end-bl+1:end,bu+1:end) = eye(bl);
else
    [m,n] = size(A);
    mmid = ceil(m/2);
    nmid = ceil(n/2);
    B.ml = mmid;
    B.nl = nmid;
    B.mr = m - mmid;
    B.nr = n - nmid;
    B.leafnode = 0;
    B.Rl = zeros(bl+bu); B.Rl(1:bl,1:bl) = eye(bl);
    B.Rr = zeros(bl+bu); B.Rr(end-bu+1:end,end-bu+1:end) = eye(bu);
    B.Wl = zeros(bl+bu); B.Wl(1:bu,1:bu) = eye(bu);
    B.Wr = zeros(bl+bu); B.Wr(end-bl+1:end,end-bl+1:end) = eye(bl);
    B.B21 = [zeros(bl,bu),full(A(mmid + 1:mmid + bl, nmid - bl + 1:nmid));zeros(bu,bl+bu)];
    B.B12 = [ zeros(bl, bl + bu);full(A(mmid-bu+1:mmid,nmid+1:nmid+bu)),zeros(bu,bl); ];
    B.A11 = hss_build_band_ric(A(1:mmid,1:nmid), bl, bu, k-1);
    B.A22 = hss_build_band_ric(A(mmid+1:end,nmid+1:end), bl, bu, k-1);
    
end
end
