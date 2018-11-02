function z = cauchy_matvec(x, y, v, tol)
% Compute the matrix vector product C(x, y) v where C_ij = 1/(xi - yj)

if ~exist('tol', 'var')
    tol = 1e-8;
end

iprec = -2;

if tol < .5
    iprec = iprec + 1;
end
if tol < .5e-1
    iprec = iprec + 1;
end
if tol < .5e-2
    iprec = iprec + 1;
end
if tol < .5e-3
    iprec = iprec + 1;
end
if tol < .5e-6
    iprec = iprec + 1;
end
if tol < .5e-9
    iprec = iprec + 1;
end
if tol < .5e-12
    iprec = iprec + 1;
end
if tol < .5e-15
    iprec = iprec + 1;
end

m = length(x);
n = length(y);

xisy = isequal(x, y);

if size(v, 1) ~= n
	error('CAUCHY_MATVEC:: dimension mismatch')
end

x = x(:).';
y = y(:).';

p = size(v, 2);
z = zeros(m, p);
for j = 1:p
	if xisy
		res = zfmm2dpart(iprec, n, [real(x); imag(x)], v(:, j), 1, 0, 0);
        z(:, j) =  res.pot.';
	else
		res = zfmm2dpart(iprec, n, [real(y); imag(y)], v(:, j), 0, 0, 0, ...
			m, [ real(x); imag(x); ], 1, 0, 0);
        	z(:, j) =  res.pottarg.';
	end
end
