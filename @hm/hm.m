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
		end
		
		if nargin == 1
			obj = create_h_matrix(varargin{1});
		end
		
		if nargin > 1
			switch varargin{1}
				case 'low-rank'
					obj = create_low_rank_h_matrix(varargin{2:end});
				otherwise
					error('Unsupported constructor mode');
			end
		end
		

		end

		
	end


end

