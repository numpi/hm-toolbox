function nrm = norm(H, nrm_type)
%NORM Estimate the norm of H

if ~exist('nrm_type', 'var')
    nrm_type = 2;
end

if ischar(nrm_type) && strcmp(nrm_type, 'fro')
    nrm = norm_frobenius(H);
    return;
end

n = size(H, 2);

nrm = normest_Afun(@(x) H * x, @(x) H' * x, n); 

end

