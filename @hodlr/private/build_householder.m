function [b, u, t] = build_householder(v)
	u = v;
	a = norm(v);
	t = a * sign(u(1));
	u(1) = u(1) + t; 	
	b = 2 / (u' * u);
	t = -t;
end
