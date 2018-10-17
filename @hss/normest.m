function nrm = normest(A, nrm_type)
if ~exist('nrm_type', 'var')
    nrm_type = 2;
end
nrm = norm(A, nrm_type);

