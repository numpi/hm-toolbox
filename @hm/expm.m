function eT = expm(A, method)
%EXPM Evaluate the matrix exponential of A. 

nrm = norm(A);
n = size(A, 2);
if nrm == 0
	eT = hm('diagonal', ones(n,1));
	return;
end
h = max(floor(log2(nrm)) + 2, 0);
A = A * (1 / 2^h);

if ~exist('method','var')
	method = 'pade';
end
if strcmp(method,'taylor')
	maxit = 12;
	eT = hm('diagonal', ones(n,1));
	
	tempT = eT;
	for i=1:maxit
		tempT = tempT * A * (1 / i);
		eT = eT + tempT;
	end
elseif strcmp(method,'pade')
	
	c = 1 / 2;
	eTn = cqt(1, 1, [], [], T.sz(1), T.sz(2)) + c*T;
	eTd = cqt(1, 1, [], [], T.sz(1), T.sz(2)) - c*T;
	
	q = 6;
	p = 1;
	X = T;
	for k = 2 : q
		c = c * (q-k+1) / (k*(2*q-k+1));
		X = T * X;
		cX = c*X;
		eTn = eTn + cX;
		if p
			eTd = eTd + cX;
		else
			eTd = eTd - cX;
		end
		p = ~p;
	end
	
	eT =  eTn * inv(eTd);
else
	error('Invalid parameter method in EXPM');
end

% eT = eT^(2^h);
for i = 1 : h
	eT = eT * eT;
end

