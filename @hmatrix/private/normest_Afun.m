function nrm = normest_Afun(Afun, Afunt, n, tol)

if ~exist('tol', 'var')
	tol = 1e-6;
end

s = 0;

v = randn(n, 1);

for i = 1 : 30
    olds = s;
	
    s = norm(v); v = v / s;
    
	% This version might be used in case we want to implement a subspace
	% iteration --- that might be more efficient for blocking. At the
	% moment the vector iteration is faster. 
	%[v, s] = qr(v, 0); s = norm(s);
	
    if abs(sqrt(olds) - sqrt(s)) < abs(s) * tol || s == 0
        break;
    end
    
    w = v;
    w = Afun(w);
    w = Afunt(w);
    v = w;
end

nrm = sqrt(s);

end
