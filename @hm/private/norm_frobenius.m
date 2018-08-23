function F = norm_frobenius(H)
%NORM_FROBENIUS Compute the Frobenius norm of H.
%
% F = NORM_FROBENIUS(H) evaluates the Frobenius norm of H. The computed
%     norm is exact, up to rounding errors.

if ~isempty(H.F)
    F = norm(H.F, 'fro');
else
    F = norm(H.A11, 'fro')^2 + norm(H.A22, 'fro')^2;
    
    % Add the contribution of the low-rank parts
    F = F + sum(sum((H.U21' * H.U21) .* (H.V21' * H.V21)));
    F = F + sum(sum((H.U12' * H.U12) .* (H.V12' * H.V12)));
    
    F = sqrt(F);
end


end

