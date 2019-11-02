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
            
            % Similarly to @hodlr and @hss, we allow the user to specify a
            % custer, which needs to be an empty @hmatrix object with just
            % the sizes and the admissibility flags set correctly.
            H = [];
            for j = 1 : length(varargin)
                if strcmp(varargin{j}, 'cluster')
                    % Build the necessary cluster, and remove the options
                    % from the input data
                    if length(varargin) < j + 1
                        error('You need to specify a cluster after the ''cluster'' keyword');
                    end
                    
                    H = varargin{j+1};
                    
                    % Filter out the cluster options from the input and
                    % proceed. 
                    varargin = varargin([ 1 : j-1, j+2:length(varargin) ]);
                    break;
                end
            end
            
            if ischar(varargin{1})
                switch varargin{1}
                    case 'handle'
                        obj = hmatrix_from_aca(H, varargin{2:4});          
                    case 'eye'
                        obj = hmatrix_eye(H, varargin{2});
                    case 'banded'
                        % At the moment we just rely on the Lanczos constructor
                        obj = hmatrix(H, varargin{2});
                    case 'low-rank'
                        obj = hmatrix_from_low_rank(H, varargin{2}, varargin{3});
                    otherwise
                        error('Unsupported constructor');
                end
            else
                % Try to construct the matrix using SVDs
                obj = hmatrix_from_full(H, varargin{1});
            end
        end                
    end
end

