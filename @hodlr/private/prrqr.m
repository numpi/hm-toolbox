function [Q, R] = prrqr(A, tol)
%

m = size(A, 1);
n = size(A, 2);

Q = eye(m);

colnrm = arrayfun(@(i) norm(A(:,i)), 1 : n);
p = 1 : n;

for j = 1 : n
	[nrm, jj] = max(colnrm(j:end));
	jj = j + jj - 1;
	
	if nrm * sqrt(n - j + 1) < tol * colnrm(1)
		j = j - 1;
		break;
	end
	
	% Swap column j and jj
	tmp = A(:,j); A(:,j) = A(:,jj); A(:,jj) = tmp;
	tmp = colnrm(j); colnrm(j) = colnrm(jj); colnrm(jj) = tmp;
	p([ j, jj ]) = p([ jj, j ]);
	
	[b, u, t] = build_householder(A(j:end,j));
	
	A(j:end,j+1:end) = A(j:end,j+1:end) - b * u * (u' * A(j:end,j+1:end));
	A(j,j) = t; A(j+1:end,j) = 0;
	Q(1:end,j:end) = Q(1:end,j:end) - b * (Q(1:end,j:end) * u) * u';
	
	colnrm = arrayfun(@(i) norm(A(j+1:end,i)), 1 : n);
end

R = A;

ip(p) = 1:n;

Q = Q(:,1:j);
R = R(1:j,ip);

end

function [b, u, t] = build_householder(v)
	u = v;
	a = norm(v);
	t = a * sign(u(1));
	u(1) = u(1) + t; 	
	b = 2 / (u' * u);
	t = -t;
end
