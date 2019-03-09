function nrm = normest(H, nrm_type)
if ~exist('nrm_type', 'var')
    nrm_type = 2;
end
nrm = norm(H, nrm_type);

