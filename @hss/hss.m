classdef hss
%HSS HSS matrices
%
% H = HSS(A) constructs an HSS representation of the matrix A, using the
%     algorithodlr described in [2]. This procedure has a cost O(n^2), where
%     n is the size of A, provided that the off-diagonal rank is negligible
%     with respect to n. 
%
%     If A is sparse, then the random sampling constructor described in
%     HSS('handle', ...) below is used.  
%
% H = HSS('banded', A) constructs an HSS representation of a banded matrix
%     A. The matrix A can be either sparse or dense. 
%
% H = HSS('banded', A, B) can be used to specify the symmetric bandwidth B 
%     of the matrix A.
%
% H = HSS('banded', A, BL, BU) specifies different lower and upper 
%     bandwidth BL and BU, respectively.
%
% H = HSS('chebfun2', F, XDOM, YDOM, M, N) constructs the M x N matrix
%     containing the samplings of the bivariate function F over a uniform
%     grid of the square XDOM x YDOM. The procedure relies on separable
%     approximation of F(X,Y) as provided by the Chebfun package.
%
% H = HSS('diagonal', D) constructs the diagonal matrix with the entries of
%     the vector D on the main diagonal. 
%
% H = HSS('eye', N) constructs an HSS representation of the N x N identity 
%     matrix. 
%
% H = HSS('handle', AFUN, AFUNT, AEVAL, M, N) constructs an HSS matrix
%     using the random sampling based algorithodlr in [1]. It requires the
%     handle function AFUN and AFUNT which perform the matrix-vector
%     products A*v and A'*v, respectively, and AEVAL which, given two
%     integer vectors I, J returns the submatrix A(I, J). M and N are the
%     number of rows and columns of A. 
%
% H = HSS('low-rank', U, V) construct an HSS representation of the low-rank
%     matrix U*V'.
%
% H = HSS('ones', M, N) constructs an HSS representation of the rank-1
%     M x N matrix of all ones. 
%
% H = HSS('toeplitz', C, R) constructs the Toeplitz matrix with C as first
%     column and R as first row. The representation is constructed using
%     the 'handle' constructor, and fast Toeplitz-vector multiplication. 
%
% H = HSS('zeros', M, N) constructs the HSS representation of the M x N
%     zero matrix. 
%
% All the constructors support an additional 'cluster' keyword that allows
% to specify custom row and column clusters. These are described as a
% vector of indices J = [J(1), ..., J(2^P)], such that the partitioning at
% the lowest level P is 
%
%        (1, J(1))    (J(1)+1, J(2))   ...   (J(2^(P-1)+1), J(2^P)), 
%
% J(2^P) = N. If J(I) = J(I+1) the corresponding leafnode is assumed to be
% missing from the tree. The cluster can be specified with the syntax
%
%   H = HSS(..., 'cluster', rowcluster, colcluster). 
%
% If colcluster is omitted then it is assumed that rowcluster == colcluster.
%
% The partitioning of an HSS matrix can be retrieved calling CLUSTER(H).  
%
%[1] Martinsson, P. G. (2011). A fast randomized algorithodlr for computing a
%    hierarchically semiseparable representation of a matrix. SIAM Journal
%    on Matrix Analysis and Applications, 32(4), 1251-1274.
%
%[2] Xia, J., Chandrasekaran, S., Gu, M., & Li, X. S. (2010). Fast 
%    algorithodlrs for hierarchically semiseparable matrices. Numerical 
%    Linear Algebra with Applications, 17(6), 953-976.
    
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
            
            rowcluster = [];
            colcluster = [];
            
            % Find the first string parameter after varargin{1}
            charpos = 2;
            while charpos <= nargin && ~ischar(varargin{charpos})
                charpos = charpos + 1;
            end
            
            if charpos <= nargin && strcmp(varargin{charpos}, 'cluster')
                rowcluster = varargin{charpos + 1};
                if nargin >= charpos + 2
                    colcluster = varargin{charpos + 2};
                else
                    colcluster = rowcluster;
                end
            end
            
            if ~ischar(varargin{1})
                A = varargin{1};
                
                if issparse(A)
                    obj = hss('handle', ...
                        @(v) A * v, @(v) A' * v, @(i,j) full(A(i,j)), ...
                        size(A, 1), size(A, 2), 'cluster', ...
                        rowcluster, colcluster);
                else
                    obj = hss_build_hss_tree(size(A, 1), size(A, 2), ...
                                hssoption('block-size'), rowcluster, ...
                                colcluster);
                    obj = hss_from_full(obj, A);
                end
                
                return;
            end
            
            if nargin > 1
                switch varargin{1}
                    case 'banded'
                        obj = hss_build_hss_tree(size(varargin{2}, 1), ...
                            size(varargin{2}, 2), hssoption('block-size'), ...
                            rowcluster, colcluster);
                        obj = hss_from_banded(obj, varargin{2:charpos-1});

                    case 'cauchy'
                        %obj = hodlr2hss(hodlr('cauchy', varargin{2:end}));
                        obj = hss_from_cauchy(varargin{2:end});

                    case 'chebfun2'
                        obj = hodlr2hss(hodlr('chebfun2', varargin{2:end}));

                    case 'diagonal'
                        obj = hss_build_hss_tree(length(varargin{2}), ...
                            length(varargin{2}), hssoption('block-size'), ...
                            rowcluster, colcluster);
                        obj = hss_build_diagonal(obj, varargin{2:charpos-1});

                    case 'eye'
                        n = varargin{2};

                        if ~check_cluster_equality(rowcluster, colcluster)
                            error('row and column cluster must match for the identity matrix');
                        end

                        obj = hss('diagonal', ones(n, 1), 'cluster', rowcluster);

                    case 'handle'
                        if charpos < 7
                            error('Unsufficient parameters for the handle constructor');
                        end
                        
                        obj = hss_build_hss_tree(varargin{5}, varargin{6}, ...
                                hssoption('block-size'), rowcluster, ...
                                colcluster);

                        obj = hss_from_random_sampling(obj, varargin{2:charpos-1});

                    case 'low-rank'
                        obj = hss_build_hss_tree(size(varargin{2}, 1), ...
                            size(varargin{3}, 1), hssoption('block-size'), ...
                            rowcluster, colcluster);
                        obj = hss_build_low_rank(obj, varargin{2:charpos-1});

                    case 'ones'
                        m = varargin{2};
                        if charpos > 3
                            n = varargin{3};
                        else
                            n = m;
                        end

                        obj = hss('low-rank', ones(m, 1), ones(n, 1), ...
                            'cluster', rowcluster, colcluster);
                        
                    case 'toeplitz'
                        if charpos == 4
                            m = length(varargin{2});
                            n = length(varargin{3});
                        elseif charpos == 5
                            m = varargin{4};
                            n = m;
                        else
                            m = varargin{4};
                            n = varargin{5};
                        end

                        obj = hss_from_symbol(varargin{2:3}, m, n, ...
                            rowcluster, colcluster);

                    case 'zeros'
                        m = varargin{2};
                        if charpos > 3
                            n = varargin{3};
                        else
                            n = m;
                        end
                        obj = hss_build_hss_tree(m, n, ...
                            hssoption('block-size'), rowcluster, colcluster);
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
