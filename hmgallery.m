function varargout = hmgallery(name, varargin)
%HMGALLERY Generate examples of H-matrices. 
%
% H = HMGALLERY('laplacian', n) generates a n x n discretization of the
% Laplacian operator on [-1, 1] using a grid of n+1 points, assuming 
% Dirichlet boundary conditions. 

switch name
    case 'laplacian'
        if isempty(varargin)
            n = 2048;
        else
            n = varargin{1};            
        end
        
        varargout{1} = hm('tridiagonal', ...
            spdiags((n - 1)^2 * ones(n,1) * [ -1 2 -1 ], -1:1, n, n));
		
case 'haber'
		if isempty('varargin')
			n = 2048;
		else
			n = varargin{1};		
		end
		
		m = 6 * n;
		
		% We first generate the matrices A and P as sparse matrices. 
		A = sparse(m, m);
		P = sparse(m, m);

		a = -1.36;
		e = 0.34;

		% We build the diagonal blocks
		for i = 1 : n
			range = ((i-1)*6+1) : i*6;
			A(range, range) = spdiags(ones(6,1) * [ e, a, e ], -1:1, 6, 6);
			P(range, range) = -.2 * ones(6) - .8 * eye(6);
		end

		% And now super and sub-diagonals
		for i = 1 : n - 1
			range = ((i-1)*6+1) : i*6;

			A(range, range + 6) = e * speye(6);
			A(range + 6, range) = e * speye(6);

			P(range, range + 6) = -.1 * ones(6);
			P(range + 6, range) = -.1 * ones(6);
		end
		
		varargout{1} = hm('banded', A, 6);
		varargout{2} = hm('banded', P, 11);
		
    otherwise
        error('Unsupported example name');
end


end

