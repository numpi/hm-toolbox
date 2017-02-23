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
        
			function [U, V] = compress_matrix(Uold)
			%COMPRESS_MATRIX Get a low-rank representation for Uold.  
			  threshold = eps;

			  [U,S,V] = svd(Uold);

			  rk = sum(diag(S) > threshold);

			  U = U(:,1:rk) * sqrt(S(1:rk,1:rk));
			  V = V(:,1:rk) * sqrt(S(1:rk,1:rk));
			end                
		end
		
		function H = create_banded_h_matrix(obj, A, band)
		%CREATE_BANDED_H_MATRIX Create a banded H-matrix. 
			H = obj;

			block_size = hmoption('block-size');

			if size(A, 1) <= block_size
				H.F = A;
				H.sz = size(A);
			else
				mp = ceil(size(A, 1) / 2);
				n = size(A, 1);

				if band <= min(n - mp)
					H.A11 = create_banded_h_matrix(hm(), A(1:mp,1:mp), band);
					H.A22 = create_banded_h_matrix(hm(), A(mp+1:end,mp+1:end), band);

					H.U12 = [ zeros(mp-band,band) ; A(mp-band+1:mp,mp+1:mp+band) ];
					H.V12 = [ eye(band) ; zeros(n - mp - band, band) ];

					H.U21 = [ A(mp+1:mp+band,mp-band+1:mp) ; zeros(n - mp - band, band) ];
					H.V21 = [ zeros(mp-band, band) ; eye(band) ];
				else
					H = create_h_matrix(H, A);
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

				H.U12 = [ zeros(mp-1,1) ; A(mp,mp+1) ];
				H.V12 = [ 1 ; zeros(n - mp - 1, 1) ];

				H.U21 = [ A(mp+1,mp) ; zeros(n - mp - 1, 1) ];
				H.V21 = [ zeros(mp-1,1) ; 1 ];

				H.sz = size(A);
			end
		end

		
	end


end

