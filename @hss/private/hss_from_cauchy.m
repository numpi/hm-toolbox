function A = hss_from_cauchy(x, y)
	if ~isvector(x) || ~isvector(y)
		error('HSS_FROM_CAUCHY:: x and y arguments must be vectors')
	end
	x = x(:); y = y(:);
	m = length(x); n = length(y);
	tol = hssoption('threshold');
	A = hss('handle', @(v) cauchy_matvec(x, y, v, tol), ...
        @(v) cauchy_matvec(-conj(y), -conj(x), v, tol), ...
        @(i, j) 1./(x(i) - y(j).'), m, n);
    %C = 1 ./ (x - y.'); C(C == inf) = 0;
    %A = hss('handle', @(v) C*v, @(v) C' * v, @(i, j) C(i,j), m, n); 
end
