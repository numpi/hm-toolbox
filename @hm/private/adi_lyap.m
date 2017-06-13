function X = adi_lyap(A, C)
%ADI_LYAP Solve the Lyapunov equation using the ADI

converged = false;

n = size(A, 1);
X = hm('diagonal', zeros(n, 1));

II = hm('diagonal', ones(n,1));

j = 0;

a = min(eig(full(A)));
b = max(eig(full(A)));

[p,q] = adi_param_syl(a, b, a, b, n);
[p,q] = adi_shift_order(logspace(log10(a), log10(b), 100), p);

while ~converged
	j = j + 1;
	pj = p(j); % Get the shift
	
	Xold = X;
	X = -2 * pj * inv(A + pj * II) * C * inv(A + pj * II)' ...
		+ inv(A + pj * II) * (A - pj * II) * X * (A - pj * II)' * inv(A + pj * II)';
	
	%j,
	%norm(Xold - X) / norm(X)
	
	if norm(Xold - X) < norm(X) * 1e-3 || j == 100
		converged = true;
	end
end

end

