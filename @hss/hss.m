classdef hss
    
    properties
        % top and bottom blocks of the matrix.
        B12
        B21
        
        % Factorization of the upper triangular block as U12 * V12'
        U
        V
        
        % Factorization of the lower triangular block as U21 * V21'
        Rl
        Rr
        Wl
        Wr
        
        % Size of the matrix
        ml
        nl
        mr
        nr
        
        % Dense version of the matrix, if the size is smaller than the
        % minimum allowed block size.
        D
        
        topnode
        leafnode
        
        A11
        A22
        
    end
    
    methods
        
        function obj = hss(varargin)
            %HSS Create a new Hierarchical matrix.
            if nargin == 0
                return;
            end
            
            obj = hss();
            
            if nargin == 1
                obj = hss_from_full(varargin{1});
                return;
            end
            
            if nargin > 1
                switch varargin{1}
                    case 'low-rank'
                        obj = hss_build_low_rank(varargin{2:end});
                    case 'diagonal'
                        obj = hss_build_diagonal(varargin{2:end});
                    case 'banded'
                        obj = hss_from_banded(varargin{2:end});
                    case 'chebfun2'
                        obj = hm2hss(hm('chebfun2', varargin{2:end}));
                    case 'cauchy'
                        obj = hm2hss(hm('cauchy', varargin{2:end}));
                        %obj = hss_from_cauchy(varargin{2:end});
                    case 'handle'
                        if length(varargin) < 6
                            error('Unsufficient parameters for the handle constructor');
                        end
                        obj = hss_from_random_sampling(varargin{2:end});
                    case 'toeplitz'
                        obj = hss_from_symbol(varargin{2:end});
                    otherwise
                        error('Unsupported constructor mode');
                end
            end
        end
        
    end
    
    %
    % Start of the private methods used to instantiate the HSS objects
    %
    methods (Access = private)
        
    end
end
