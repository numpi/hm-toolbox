function [p,q] = adi_param_syl(a,b,c,d,n)
%ADI_PARAM_SYL
%Computes the optimal parameters pj, qj for the general
%ADI minimax problem
%min_{p1,...pr; q1,...qr} ...
%max{x in [a,b], y in [c,d]} ...
%prod_{j=1}^r abs( (x-qj)(y-pj) / (x+pj)(y+qj) )
%These parameters are used in the ADI method for solving the Sylvester
%equation A X + X B + C D’ = 0. The parameters provided by this function
%are optimal if A=A’ has eigenvalues in [a,b] and B=B’ has eigenvalues in
%[c,d].
%INPUT:
%a,b,c,d - intervals [a,b] and [c,d]
%n - number of parameters desired, an integer >= 1
%OUTPUT:
%p,q - an n by 1 vectors of the parameters
sgn = sign(a);
if (((b/a) > 1e8) || ((d/c) > 1e8))
	p = logspace(log10(a),log10(b),n);
	q = logspace(log10(c),log10(d),n);
else
	mh = 2*(b-a)*(d-c)/(a+c)/(b+d);
	kp = 1/(1 + mh + sqrt(mh*(mh+2)));
	if (a==c) && (b==d)
		g = 0;
	else
		g = 2*(kp*(b+d)-(a+c))/((a+c)*(b-d)+kp*(b+d)*(c-a));
	end
	f = (2+g*(b-d))/(b+d);
	h = kp*(c-a+2*a*c*g)/(a+c);
	m = 1-kp^2;
	K = ellipke(m);
	u = (2*(1:n)'-1)/(2*n) * K;
	[~,~,dn] = ellipj(u,m);
	w = dn;
	p = (h+w)./(w*g+f);
	q = (h-w)./(w*g-f);
end
p = sgn*p;
q = sgn*q;