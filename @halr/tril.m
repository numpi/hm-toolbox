function tH = tril(H)
%TRIL computes the lower triangular part of a HALR matrix

tH = H;

if is_leafnode(H)
    if H.admissible
        % In this case it is very likely that the block will end up being
        % unstructured
        tH.admissible = false;
        tH.F = tril(H.U * H.V');
    else
        tH.F = tril(H.F);
    end
else
    tH.A21 = H.A21;
    tH.A11 = tril(H.A11);
    tH.A22 = tril(H.A22);
    
    % Zero the upper off-diagonal block
    [m1, ~] = size(H.A11);
    [~, n2] = size(H.A22);
    
    % Zero off-diagonal blocks
    tH.A12 = halr;
    tH.A12.sz = [m1 n2];
    tH.A12.U = zeros(m1, 0);
    tH.A12.V = zeros(n2, 0);
    tH.A12.admissible = true;
end
