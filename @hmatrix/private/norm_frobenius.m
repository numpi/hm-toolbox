function F = norm_frobenius(H)
%NORM_FROBENIUS Compute the Frobenius norm of H.
%
% F = NORM_FROBENIUS(H) evaluates the Frobenius norm of H. The computed
%     norm is exact, up to rounding errors.

if is_leafnode(H)
    if H.admissible
        F = sqrt(sum(sum((H.U' * H.U) .* (H.V' * H.V))));
    else
        F = norm(H.F, 'fro');
    end
else
    F = norm(H.A11, 'fro')^2 + norm(H.A22, 'fro')^2 + ...
        norm(H.A12, 'fro')^2 + norm(H.A21, 'fro')^2;
    F = sqrt(F);
end


end

