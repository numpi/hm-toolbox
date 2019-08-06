function [HL, HU] = hodlr_lu(H)
%HODLRATRIX_LU LU factorization for hodlr objects

HL = H;
HU = H;

if size(H, 1) ~= size(H, 2)
    error('LU is supported only for square HODLR matrices');
end

HL.sz = [size(H, 1), size(H, 1)];

if is_leafnode(H)
    [HL.F, HU.F] = lu(H.F);
else
    m1 = H.A11.sz(1);
    n1 = H.A11.sz(2);
    
    n = H.sz(2);
    m = H.sz(1);
    
    HL.U12 = zeros(m1,0);
    HL.V12 = zeros(m-m1,0);
    HU.U21 = zeros(m-m1,0);
    HU.V21 = zeros(n1,0);
    
    [HL.A11, HU.A11] = hodlr_lu(H.A11);
    HU.U12 = solve_lower_triangular(HL.A11, H.U12);
    % HL.V21 = (H.V21' / HU.A11)';
    HL.V21 = solve_upper_triangular2(HU.A11, H.V21')';
    [HL.A22, HU.A22] = hodlr_lu(hodlr_rank_update(H.A22,-HL.U21*(HL.V21'*HU.U12),HU.V12));
end


end

