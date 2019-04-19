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
v = randn(n, 1);

Ht = H';

s = 0;

for i = 1 : 10
    olds = s;
    s = norm(v);
    
    if abs(olds - s) < abs(s) * 1e-3 || s == 0
        break;
    end
    
    v = v / s;
    w = v;
    w = H * w;
    w = Ht * w;
    v = w;
end

nrm = sqrt(s);

end

