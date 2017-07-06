function X = lyap(A, C, varargin)
%LYAP Solve the Lyapounov equation AX + XA' + C = 0

N = 32;

[x, w] = hm_legpts(N);

p = inputParser;

addOptional(p, 'debug', false, @islogical);
addOptional(p, 'expm',  'pade', @ischar);
addOptional(p, 'method', 'expm', @ischar);

parse(p, varargin{:});

debug = p.Results.debug;
expm_method = p.Results.expm;
method = p.Results.method;

if strcmp(method, 'sign')
	X = sign_lyap(A, C);
	return;
end

if strcmp(method, 'adi')
	X = adi_lyap(A, C, 1e-6, 100);
	return;
end

% Acceleration parameter. 
L = 5;

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

% X = partialSum(xx, ww, @f);
X = ww(1) * f(xx(1));
for i = 2 : length(ww)
    X = X + ww(i) * f(xx(i));
end

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
      if strcmp(expm_method, 'mixed')
          if abs(x) * nrm > 8
              expm_method = 'ratcheb';
          else
              expm_method = 'pade';
          end
      end
          
      eA = expm(-x*A, expm_method, 6, nrm * abs(x));
	  nevals = nevals + 1;
	  Y = -eA * C * eA';
	  
	  if debug
		fprintf (' :: %d out of %d evaluations performed\n', nevals, N);
	  end
	end

end

