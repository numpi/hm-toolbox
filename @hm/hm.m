classdef hm
	
	properties
		% top and bottom blocks of the matrix. 
		A11
		A22
		
		% Factorization of the upper triangular block as U12 * V12'
		U12
		V12
		
		% Factorization of the lower triangular block as U21 * V21'
		U21
		V21
		
		% Size of the matrix
		sz
		
		% Dense version of the matrix, if the size is smaller than the
		% minimum allowed block size. 
		F
		
		% If this flag is true no compression is performed, until
		% hm_recompress is called. 
		compression_disabled
	end
	
	methods
		
		function obj = hm(varargin)
		%HM Create a new Hierarchical matrix. 		
			if nargin == 0
				obj.F = [];
				obj.sz = [0 0];
				return;
			end

			if nargin == 1
				obj = create_h_matrix(obj, varargin{1});
				return;
			end

			if nargin > 1    
				switch varargin{1}
					case 'low-rank'
						obj = create_low_rank_h_matrix(obj, varargin{2:end});
					case 'tridiagonal'
						obj = create_tridiagonal_h_matrix(obj, varargin{2});
					case 'banded'
						obj = create_banded_h_matrix(obj, varargin{2:end});
					case 'diagonal'
						obj = create_diagonal_h_matrix(obj, varargin{2:end});
					case 'sparse'
						obj = create_sparse_h_matrix(obj, varargin{2});
                    case 'chebfun2'
                        obj = create_chebfun2_h_matrix(obj, varargin{2:end});
					case 'toeplitz'
						obj = create_toeplitz_h_matrix(obj, varargin{2:end});
                    case 'cauchy'
                        warning('The CAUCHY constructor is not (yet) efficiently implemented');
                        obj = create_cauchy_h_matrix(obj, varargin{2:end});
					otherwise
						error('Unsupported constructor mode');
				end
			end
		end
		
	end
	
	% 
	% Start of the private methods used to instantiate the HM objects
	%
	
	methods (Access = private)
		
		function obj = create_h_matrix(obj, A)
			%CREATE_H_MATRIX Given a dense matrix A, construct a hierarchical
			%representation for it. 
			min_block_size = hmoption('block-size');

			obj.F = [];
			obj.sz = size(A);

			if size(A, 1) <= min_block_size && size(A, 2) <= min_block_size
				obj.F = A;
			else
				% Get the middle point
				mp = ceil(size(A, 1) / 2);

				obj.A11 = create_h_matrix(hm(), A(1:mp,1:mp));
				obj.A22 = create_h_matrix(hm(), A(mp+1:end,mp+1:end));

				[obj.U21, obj.V21] = compress_matrix(A(mp+1:end,1:mp));
				[obj.U12, obj.V12] = compress_matrix(A(1:mp,mp+1:end));
            end                       
		end
		
		function obj = create_sparse_h_matrix(obj, A)
		%CREATE_SPARSE_H_MATRIX Create an H-matrix starting from a sparse one
		
			n = size(A, 2);
			block_size = hmoption('block-size');
			
			if n <= block_size
				obj.F = A;
			else
				mp = ceil(size(A, 1));
				obj.A11 = create_sparse_h_matrix(hm(), A(1:mp, 1:mp));
				obj.A22 = create_sparse_h_matrix(hm(), A(mp+1:end,mp+1:end));
				
				% FIXME: Missing implementation of two-sided Lanczos method
			end
		end
		
		function H = create_banded_h_matrix(obj, A, band)
		%CREATE_BANDED_H_MATRIX Create a banded H-matrix. 
			H = obj;

			block_size = hmoption('block-size');

			if size(A, 1) <= block_size
				H.F = full(A);
				H.sz = size(A);
			else
				mp = ceil(size(A, 1) / 2);
				n = size(A, 1);

				if band <= min(n - mp)
					H.A11 = create_banded_h_matrix(hm(), A(1:mp,1:mp), band);
					H.A22 = create_banded_h_matrix(hm(), A(mp+1:end,mp+1:end), band);

					H.U12 = [ zeros(mp-band,band) ; full(A(mp-band+1:mp,mp+1:mp+band)) ];
					H.V12 = [ eye(band) ; zeros(n - mp - band, band) ];

					H.U21 = [ full(A(mp+1:mp+band,mp-band+1:mp)) ; zeros(n - mp - band, band) ];
					H.V21 = [ zeros(mp-band, band) ; eye(band) ];
                    
                    % Perform a compression
                    [H.U21, H.V21] = compress_factors(H.U21, H.V21, norm(H.U21, 'fro'));                    
                    [H.U12, H.V12] = compress_factors(H.U12, H.V12, norm(H.U12, 'fro'));
				else
					H = create_h_matrix(H, full(A));
				end

				H.sz = size(A);
			end
		end
		
		function obj = create_diagonal_h_matrix(obj, D)
		%CREATE_DIAGONAL_H_MATRIX Create an H-matrix with the specified diagonal. 
			n = length(D);

			if length(D) <= hmoption('block-size')
				obj.F = diag(D);
				obj.sz = [ n, n ];
			else
				mp = ceil(length(D) / 2);
				obj.A11 = create_diagonal_h_matrix(hm(), D(1:mp));
				obj.A22 = create_diagonal_h_matrix(hm(), D(mp+1:end));
				obj.U12 = zeros(mp, 0);
				obj.V12 = zeros(n - mp, 0);
				obj.U21 = zeros(n - mp, 0);
				obj.V21 = zeros(mp, 0);

				obj.sz = [ n, n ];
			end
		end
		
		function obj = create_low_rank_h_matrix(obj, U, V)
		%CREATE_LOW_RANK_H_MATRIX Create a low rank H matrix. 
			global hm_block_size

			if isempty(hm_block_size)
				hmoption('block-size');
			end

			obj.sz = [ size(U, 1), size(V, 1) ];

			if obj.sz(1) <= hm_block_size
				obj.F = U * V';
			else
				mp = ceil(obj.sz(1) / 2);
				obj.A11 = create_low_rank_h_matrix(hm(), U(1:mp,:), V(1:mp,:));
				obj.A22 = create_low_rank_h_matrix(hm(), U(mp+1:end,:), V(mp+1:end,:));
				obj.U12 = U(1:mp,:);
				obj.V12 = V(mp+1:end,:);
				obj.U21 = U(mp+1:end,:);
				obj.V21 = V(1:mp,:);
			end
		end		
		
		function H = create_tridiagonal_h_matrix(obj, A)
		%CREATE_TRIDIAGONAL_H_MATRIX Create a tridiagonal H-matrix

			H = obj;

			block_size = hmoption('block-size');

			if size(A, 1) <= block_size
				H.F = full(A);
				H.sz = size(A);
			else
				mp = ceil(size(A, 1) / 2);
				n = size(A, 1);

				H.A11 = create_tridiagonal_h_matrix(hm(), A(1:mp,1:mp));
				H.A22 = create_tridiagonal_h_matrix(hm(), A(mp+1:end,mp+1:end));

				H.U12 = [ zeros(mp-1,1) ; full(A(mp,mp+1)) ];
				H.V12 = [ 1 ; zeros(n - mp - 1, 1) ];

				H.U21 = [ full(A(mp+1,mp)) ; zeros(n - mp - 1, 1) ];
				H.V21 = [ zeros(mp-1,1) ; 1 ];

				H.sz = size(A);
			end
        end
        
        function obj = create_chebfun2_h_matrix(obj, fct, xdom, ydom, n)
            block_size = hmoption('block-size');
            obj.sz = [ n, n ];

            if ~exist('chebfun2', 'class')
                error('Chebfun not found: did you forget to add it to the path?');
            end
            
            x = linspace(xdom(1), xdom(2), n);
            y = linspace(ydom(1), ydom(2), n);
            
            if n <= block_size
                obj.F = fct( ones(n,1) * x, y.' * ones(1,n) );
            else                
				mp = ceil(n / 2);
                obj.A11 = create_chebfun2_h_matrix(obj(), fct, ...
                    [ x(1), x(mp) ], ...
                    [ y(1), y(mp) ], mp);
                obj.A22 = create_chebfun2_h_matrix(hm(), fct, ...
                    [ x(mp+1), x(end) ], ...
                    [ y(mp+1), y(end) ], n - mp);
                
                % Create the low-rank block A12 and A21
                [obj.U12, obj.V12] = chebfun2_low_rank(fct, ...
                    [ x(mp+1), x(end) ], ...
                    [ y(1), y(mp) ], ...
                    mp, n - mp);
                [obj.U21, obj.V21] = chebfun2_low_rank(fct, ...
                    [ x(1), x(mp) ], ...
                    [ y(mp+1), y(end) ], ...
                    n - mp, mp);
			end
			
			obj = compress_hmatrix(obj);
		end
		
		function obj = create_toeplitz_h_matrix(obj, am, ap, n)
			block_size = hmoption('block-size');
			obj.sz = [ n, n ];
			
			if length(am) > n
				am = am(1:n);
				warning('Negative symbol is longer than n: truncation has been applied');
			end
			
			if length(ap) > n
				ap = ap(1:n);
				warning('Positive symbol is longer than n: truncation has been applied');
			end			
			
			if length(am) < n
				am(n) = 0;
			end
			
			if length(ap) < n
				ap(n) = 0;
			end			
			
			if n <= block_size
				obj.F = toeplitz(am, ap);
			else
				% obj.A11 = create_toeplitz_h_matrix(hm(), am, ap, block_size);
				% obj.A22 = obj.A11;
				
				mp = ceil(n / 2);
				np = mp;
								
				tol = hmoption('threshold');
				
				% Here the truncation is made relatively to the norm of the
				% off-diagonal blocks: we might want to do it relatively to
				% the norm of the big matrix. 
				[tU21,S21,tV21] = lanczos_svd(@(v,trasp) toepmult_afun(...
					am(mp+1:end), am(mp+1:-1:2), ...
					n - mp, np, v, trasp), n - mp, np, tol);
				[tU12,S12,tV12] = lanczos_svd(@(v,trasp) toepmult_afun(...
					ap(np+1:-1:2), ap(np+1:end), ...
					mp, n - np, v, trasp), mp, n - np, tol);
				
				tU21 = tU21 * sqrt(S21);
				tV21 = tV21 * sqrt(S21);
				tU12 = tU12 * sqrt(S12);
				tV12 = tV12 * sqrt(S12);
				
				obj = initialize_toeplitz_h_matrix(obj, am, ap, n, ...
					tU12, tV12, tU21, tV21);
			end
			
		end
		
		function obj = initialize_toeplitz_h_matrix(obj, am, ap, n, tU12, tV12, tU21, tV21)
			obj.sz = [ n n ];
			
			if n <= hmoption('block-size')
				obj.F = toeplitz(am(1:n), ap(1:n));
			else
				mp = ceil(n / 2);
				np = ceil(n / 2);
				
				obj.U21 = tU21(1:(n-mp), :);
				obj.V21 = tV21(end-np+1:end,:);
				obj.U12 = tU12(end-mp+1:end,:);
				obj.V12 = tV12(1:(n-np), :);
				
				obj.A11 = initialize_toeplitz_h_matrix(hm(), am, ap, mp, tU12, tV12, tU21, tV21);				
				obj.A22 = initialize_toeplitz_h_matrix(hm(), am, ap, n-mp, tU12, tV12, tU21, tV21);
			end
        end
        
        function obj = create_cauchy_h_matrix(obj, x, y)
            n = length(x);
            
            obj.sz = [ n n ];
            
            x = reshape(x, length(x), 1);
            y = reshape(y, length(y), 1);
            
            if n <= hmoption('block-size')
                obj.F = 1 ./ (x + y.');
                for i = 1 : n
                    obj.F(i,i) = 0;
                end
            else
                mp = ceil(n / 2);
				np = ceil(n / 2);
                
                [obj.U21, obj.V21] = cauchy_lr(x(mp+1:end), y(1:np));
                [obj.U12, obj.V12] = cauchy_lr(x(1:mp), y(np+1:end));
                
                obj.A11 = hm();
                obj.A22 = hm();
                
                obj.A11 = create_cauchy_h_matrix(obj.A11, x(1:mp), y(1:np));
                obj.A22 = create_cauchy_h_matrix(obj.A22, x(mp+1:end), y(np+1:end));
            end
        end

		
    end

end

