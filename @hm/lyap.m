function X = lyap(A, C)
%LYAP Solve the Lyapounov equation AX + XA' = C

n = size(A, 1);
[x, w] = hm_legpts(32);

X = hm('diagonal',zeros(n,1));
L = 100;

nrm = norm(A);

for i=1:length(x)
	% cot(x(i)/2)^2*L
	X = X + 2*L*w(i) * f(L*cot(x(i)/2)^2) *(sin(x(i))/(1-cos(x(i)))^2);
end

	function Y = f(x)	  
	  eA = expm(-x*A, 'pade', 12, nrm * abs(x));
	  Y = eA * C * eA';	  
	end

end

