function [nrm, A] = norm_frobenius(A)
A = hss_proper(A);
nrm = sqrt(norm_frobenius_ric(A));
end

function nrm = norm_frobenius_ric(A)
if A.leafnode == 1
    nrm = norm(A.D,'fro')^2;
else
    nrm = sum(sum(A.B12.^2)) + sum(sum(A.B21.^2)) + norm_frobenius_ric(A.A11) + norm_frobenius_ric(A.A22);
end
end
