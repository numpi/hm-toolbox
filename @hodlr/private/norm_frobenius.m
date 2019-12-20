function F = norm_frobenius(H)
%NORM_FROBENIUS Compute the Frobenius norm of H.
%
% F = NORM_FROBENIUS(H) evaluates the Frobenius norm of H. The computed
%     norm is exact, up to rounding errors.

if is_leafnode(H)
    F = norm(H.F, 'fro');
else
    F11 = norm(H.A11, 'fro');
    F22 = norm(H.A22, 'fro');
    
    % Contribution from the block (2,1)
    [~, RU] = qr(H.U21, 0); %[~, RV] = qr(H.V21, 0);
    %F21 = norm(RU * RV', 'fro');
    F21 = norm(H.V21 * RU, 'fro');
    
    % Contribution from the block (1,2)
    [~, RU] = qr(H.U12, 0); %[~, RV] = qr(H.V12, 0);
    %F12 = norm(RU * RV', 'fro');
    F12 = norm(H.V12 * RU, 'fro');
    
    F = norm([ F11, F12, F21, F22 ]);
end


end

