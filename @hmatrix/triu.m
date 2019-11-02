function tH = triu(H)
%TRIL computes the upper triangular part of a HMATRIX matrix

tH = H;

if is_leafnode(H)
    if H.admissible
        % In this case it is very likely that the block will end up being
        % unstructured
        tH.admissible = false;
        tH.F = triu(H.U * H.V');
    else
        tH.F = triu(H.F);
    end
else
    tH.A12 = H.A12;
    tH.A11 = triu(H.A11);
    tH.A22 = triu(H.A22);
    
    % Zero the upper off-diagonal block
    [~, n1] = size(H.A11);
    [m2, ~] = size(H.A22);
    
    % Zero off-diagonal blocks
    tH.A21 = hmatrix;
    tH.A21.sz = [m2 n1];
    tH.A21.U = zeros(m2, 0);
    tH.A21.V = zeros(n1, 0);
    tH.A21.admissible = true;
end
