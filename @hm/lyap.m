function X = lyap(A, C, varargin)
%LYAP Solve the Lyapounov equation AX + XA' + C = 0

N = 32;

[x, w] = hm_legpts(N);

p = inputParser;

addOptional(p, 'debug', false, @isboolean);
addOptional(p, 'expm',  'pade', @ischar);

parse(p, varargin{:});

debug = p.Results.debug;
expm_method = p.Results.expm;

% Acceleration parameter. 
L = 100;

% The computation of the norm is required only for scaling and squaring
% methods, therefore we can save some time in the ratcheb case. 
if ~strcmp(expm_method, 'ratcheb')
    nrm = norm(A);
else
    nrm = 0.0;
end

% Perform the change of variables and weights
ww = (2 * L * w) .* ( sin(x) ./ (1 - cos(x)).^2 );
xx = L * cot(x ./ 2).^2;

nevals = 0;

X = partialSum(xx, ww, @f);

	function X = partialSum(xx, ww, ef)
		%PARTIALSUM Perform the sum of the function EF at some points. 
		%
		% X = PARTIALSUM(XX, WW, EF) evaluates EF at the points in XX, and
		% combines the evaluations using the weights in WW. 
		if length(xx) == 1
			X = ww * ef(xx);
		else
			mp = ceil(length(xx) / 2);
			
			X = partialSum(xx(1:mp), ww(1:mp), ef) + ...
				partialSum(xx(mp+1:end), ww(mp+1:end), ef);
		end			
	end

	function Y = f(x)	  
	  eA = expm(-x*A, expm_method, 8, nrm * abs(x));
	  nevals = nevals + 1;
	  Y = -eA * C * eA';
	  
	  if debug
		fprintf (' :: %d out of %d evaluations performed\n', nevals, N);
	  end
	end

end

