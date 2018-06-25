function X = lyap(A, C, varargin)
%LYAP Solve the Lyapounov equation AX + XA' + C = 0
%
% X = LYAP(A, C) solves the Lyapunov equation using a divide and conquer
%     method. 
%
% X = LYAP(A, C, SA) uses the D&C method exploiting the sparse version of
%     A, passed as the argument SA, to construct the Krylov spaces. 
%
% X = LYAP(A, C, 'method', METHOD) uses the method specified. The available
%     options are: 
%      - 'd&c': The D&C approach of [2]. 
%      - 'sign': The sign function iteration, see [1,3] (works only if
%         A has eigenvalue on the open left or right half-planes).
%      - 'expm': Uses a quadrature formula to approimate the solution [3]. 
%         Works for positive definite matrices. 
% 
% [1] Grasedyck, L., Hackbusch, W., & Khoromskij, B. N. (2003). Solution
%     of large scale algebraic matrix Riccati equations by use of 
%     hierarchical matrices. Computing, 70(2), 121-165.
% 
% [2] Kressner, D., Massei, S., & Robol, L. (2017). Low-rank updates 
%     and a divide-and-conquer method for linear matrix equations. 
%     arXiv preprint arXiv:1712.04349.
%
% [3] Massei, S., Palitta, D., & Robol, L. (2017). Solving rank structured
%     Sylvester and Lyapunov equations. arXiv preprint arXiv:1711.05493.

N = 64;

[x, w] = hm_legpts(N);

if ~isempty(varargin) && ~ischar(varargin{1})
	nrmA = norm(A);
	
	X = sparse_dac_lyap(A / nrmA, A' / nrmA, C / nrmA , ...
		varargin{1} / nrmA, varargin{2} / nrmA);
	return;
end

p = inputParser;

addParameter(p, 'debug', false, @islogical);
addParameter(p, 'expm',  'pade', @ischar);
addParameter(p, 'method', 'd&c', @ischar);
addParameter(p, 'parallel', false, @islogical);

parse(p, varargin{:});

debug = p.Results.debug;
expm_method = p.Results.expm;
method = p.Results.method;
parallel = p.Results.parallel;

if strcmp(method, 'sign')
	X = sign_lyap(A, C);
	return;
end

if strcmp(method, 'd&c')
	X = dac_lyap(A,A',C);
    X = compress_hmatrix(X);
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
	
	if nrm > 5.37	
		A = A / (nrm / 5.37);
		C = C / (nrm / 5.37);
		nrm = 5.37;
	end
else
    nrm = 0.0;
end

% Perform the change of variables and weights
ww = (2 * L * w) .* ( sin(x) ./ (1 - cos(x)).^2 );
xx = L * cot(x ./ 2).^2;

nevals = 0;

% X = partialSum(xx, ww, @f);
if ~parallel
	X = ww(1) * f(xx(1));
	for i = 2 : length(ww)
		X = X + ww(i) * f(xx(i));
	end
else
	n = size(A, 1);
	if exist('gcp', 'file')
		pool = gcp;
		nworkers = pool.NumWorkers * 6;		
	else
		nworkers = matlabpool('size');
	end
	
	X = hm('diagonal', zeros(n, 1));
	
	XX = cell(1, nworkers);
	for j = 1 : ceil(N / nworkers)
		base = (j - 1) * nworkers + 1;
		
		parfor i = 1 : min(length(ww) - base, nworkers)
			if strcmp(expm_method, 'mixed')		  
				if abs(xx(i+base)) * nrm > 64
					expm_l_method = 'ratcheb';
				else
					expm_l_method = 'pade';
				end
			else
				expm_l_method = expm_method;
			end

			eA = expm(-xx(i+base)*A, expm_l_method, 13, nrm * abs(xx(i+base)));
			nevals = nevals + 1;
			XX{i} = -ww(i+base) * eA * C * eA';
		end
		
		for i = 1 : nworkers
			X = X + XX{i};
		end		
	end
	
	
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
          if abs(x) * nrm > 64
              expm_l_method = 'ratcheb';
          else
              expm_l_method = 'pade';
		  end
	  else
		  expm_l_method = expm_method;
	  end	  	  
	  
	  %if strcmp(expm_l_method, 'ratcheb')
	  %	  fprintf('r ');
	  %else
	  %	  fprintf('p ');
	  %end
          
      eA = expm(-x*A, expm_l_method, 13, nrm * abs(x));
	  nevals = nevals + 1;
	  Y = -eA * C * eA';
	  
	  if debug
		fprintf (' :: %d out of %d evaluations performed\n', nevals, N);
	  end
	end

end

