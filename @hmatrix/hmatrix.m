classdef hmatrix
    %HMATRIX Hierarchical matrix with arbitrary partition
    %
    % 
    
    properties
        U
        V
        F
        A11
        A12
        A21
        A22
	admissible
	sz
    end
    
    methods
        function obj = hmatrix(varargin)
            obj.admissible = false;
            
            if nargin == 0
                return;
            end            
            
            if ischar(varargin{1})
                switch varargin{1}
                    case 'handle'
                        obj = hmatrix_from_aca(varargin{2:4});          
                    case 'eye'
                        obj = hmatrix_eye(varargin{2});
                    case 'banded'
                        % At the moment we just rely on the Lanczos constructor
                        obj = hmatrix(varargin{2});
                    case 'low-rank'
                        obj = hmatrix_from_low_rank(varargin{2}, varargin{3});
                    otherwise
                        error('Unsupported constructor');
                end
            else
                % Try to construct the matrix using SVDs
                obj = hmatrix_from_full(varargin{1});
            end
        end                
    end
end

