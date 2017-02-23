function X = lyap(A, C)
%LYAP Solve the Lyapounov equation AX + XA' = C

n = size(A, 1);
[x, w] = legpts(64,[0,pi]);

%hh = ceil(log2(1e5 * norm(A))) + 1;
%AA = cell(1, 12);
%AA{1} = -A * (1 / 2^hh);
%for i = 2 : 12
%	AA{i} = AA{i-1} * AA{1} * (1 / i);
%end

X = hm('diagonal',zeros(n,1));
L = 100;
for i=1:length(x)
	% x(i)
	X = X + 2*L*w(i) * f(L*cot(x(i)/2)^2) *(sin(x(i))/(1-cos(x(i)))^2) ; 
end

	function Y = f(x)
	  if (abs(x) < 1e5) && false
		 eA = smartexpm(AA, x);
	  else		  
		eA = expm(-x*A,'taylor');
	  end
	  Y = eA *C * eA';
	end

	function Y = smartexpm(AA, t)		
		Y = hm('diagonal', ones(n, 1));
		for ii = 1 : length(AA)
			Y = Y + t^ii * AA{ii};
		end
		for ii = 1 : hh
			Y = Y * Y;
		end
	end

end

