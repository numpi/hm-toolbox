function [p,q,shiftorder,maxr] = adi_shift_order(z1,p,z2,q)
if nargin < 3
	z2 = -z1;
	q = -p;
end
n = length(p);
shiftorder = zeros(n,1);

R = zeros(length(z1),n);
for j = 1:n
	R(:,j) = abs( ( z2-q(j) ).*(z1-p(j))./(z2 + p(j)')./(z1+q(j)'));
end

for k=n:-1:2
	shift_indices = setdiff(1:n,shiftorder(k+1:end));
	maxr = zeros(k,1);
	for j=1:k
		use_these_shifts = shift_indices;
		use_these_shifts(j) = [];
		r = prod(R(:,use_these_shifts),2);
		maxr(j) = max(r);
	end
	[minmaxr,ii] = min(maxr);
	shiftorder(k) = shift_indices(ii);
end
shiftorder(1) = setdiff(1:n,shiftorder(2:end));
p = p(shiftorder);
q = q(shiftorder);
%[r, maxr] = adi_eval_complex(z1,p);

