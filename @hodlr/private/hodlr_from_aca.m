function H = hodlr_from_aca(H, Aeval, m1, m2, n1, n2)
%HODLR_FROM_ACA Constructs the HODLR blocks using Adaptive Cross

if is_leafnode(H)
    H.F = Aeval((m1:m2)', n1:n2);
else
    mm = m1 + H.A11.sz(1);
    nm = n1 + H.A11.sz(2);
    
    H.A11 = hodlr_from_aca(H.A11, Aeval, m1, mm - 1, n1, nm - 1);
    H.A22 = hodlr_from_aca(H.A22, Aeval, ...
            mm, m2, nm, n2);
        
    % ACA for the block (1,2)
    [H.U12, H.V12] = aca(@(i,j) Aeval(m1 + i - 1, nm + j - 1), size(H.A11, 1), ...
            size(H.A22,2), hodlroption('threshold'), 1);
    % ACA for the block (2,1)
    [H.U21, H.V21] = aca(@(i,j) Aeval(mm + i - 1, n1 + j - 1), size(H.A22, 1), ...
            size(H.A11, 2), hodlroption('threshold'), 1);
end

end

