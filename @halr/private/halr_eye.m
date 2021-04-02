function H = halr_eye(H, n)
%HALR_EYE Identity matrix

if isempty(H)
    H = halr_build_default_tree(n, n, halroption('block-size'));
end

H = halr_eye_ric(H);

end

function H = halr_eye_ric(H)
if is_leafnode(H)
    H.admissible = false;
    H.F = eye(size(H));
else
    H.A11 = halr_eye_ric(H.A11);
    H.A22 = halr_eye_ric(H.A22);
    
    [m1, n1] = size(H.A11);
    [m2, n2] = size(H.A22);
    
    % Zero off-diagonal blocks
    H.A12 = halr;
    H.A12.sz = [m1 n2];
    H.A12.U = zeros(m1, 0);
    H.A12.V = zeros(n2, 0);
    H.A12.admissible = true;
    
    H.A21 = halr;
    H.A21.sz = [m2 n1];
    H.A21.U = zeros(m2, 0);
    H.A21.V = zeros(n1, 0);
    H.A21.admissible = true;
end
end