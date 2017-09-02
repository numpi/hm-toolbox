function nrm = norm(A, nrm_type)
%NORM Estimate the norm of A

if ~exist('nrm_type', 'var')
    nrm_type = 2;
end

% Not yet implemented
%if ischar(nrm_type) && strcmp(nrm_type, 'fro')
%    nrm = norm_frobenius(A);
%   return;
%nd

n = size(A, 2);
v = randn(n, 1);

At = A';

s = 0;

for i = 1 : 30
	olds = s;
	s = norm(v);
	
	if abs(olds - s) < abs(s) * 1e-3
		break;
	end
	
	v = v / s;
	w = v;
	w = A * w;
	w = At * w;
	v = w;
end

nrm = sqrt(s);

end

