classdef hmatrix
    %HMATRIX Hierarchical matrix with arbitrary partition
    %
    % This class represents Hierarchical matrices in H-format, as
    % implemented by the HODLR class, but with an arbitrary (binary)
    % partitioning. 
    %
    % If no particular choice of partitioning is given, the constructor
    % chooses a default one which has a dense tridiagonal structure with
    % low-rank blocks that get larger and larger away from the main block
    % diagonals. 
    %
    % H = HMATRIX(A) constructs an HMATRIX representation of the matrix A,
    %   which might be dense or sparse, with the default partitioning. 
    %
    % To specify an arbitrary partitioning, the user has to construct a
    % quad-tree using the HMATRIX structure itself, in the following way:
    %
    %  cluster = hmatrix;
    %  cluster.sz [ m n ];
    %  cluster.admissible = true / false;
    %  
    % which sets the size and the admissibility (i.e., if it is low-rank)
    % of a certain block. Then, it is possible to attach recursively other
    % children to it to represent lower levels, using the fields A11, A12,
    % A21, and A22. For instance, a 2-level partitioning with dense
    % diagonal blocks can be achieved as follows:
    %
    %  c = hmatrix; c.sz = [ m n ];
    %  c.A11 = hmatrix; c.A11.sz = [ m1 n1 ]; c.A11.admissible = false;
    %  c.A22 = hmatrix; c.A22.sz = [ m-m1, n-n1 ]; c.A22.admissible = false;
    %  c.A21 = hmatrix; c.A22.sz = [ m-m1, n1 ]; c.A22.admissible = true;
    %  c.A12 = hmatrix; c.A12.sz = [ m1, n-n1 ]; c.A12.admissible = true;
    % 
    % To construct an HMATRIX with prescribed structure, one can then use
    % the command:
    %
    %  H = HMATRIX(A, 'cluster', c);
    %
    % The same syntax can be used with all the other constructors, listed
    % here:
    %
    % H = HMATRIX('banded', A) constructs a representation of a banded
    %     matrix given in sparse format. 
    %
    % H = HMATRIX('handle', Afun, m, n) builds an M x N matrix by evaluting
    %     the handle function Afun(I, J) at some pivots by means of ACA
    %     (Adaptive Cross Approximation). 
    
    properties
        % Left factor of the low-rank factorization of this block, in case
        % the admissible flag is set to true, and this is a leaf node. 
        U
                
        % Right factor of the low-rank factorization of this block, in case
        % the admissible flag is set to true, and this is a leaf node
        V
        
        % Dense representation of this block, in case this is a
        % non-admissible leaf node
        F
        
        % Children in position (1,1), if this is not a leaf node
        A11
        
        % Children in position (1,2), if this is not a leaf node
        A12
        
        % Children in position (2,1), if this is not a leaf node
        A21
        
        % Children in position (2,2), if this is not a leaf node
        A22
        
        % This flag if set to true on leaf nodes if they can be
        % approximated by low-rank blocks. 
        admissible
        
        % This is an 1 x 2 matrix containing the dimension of H, that is we
        % have H.sz = [m n]. 
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
                % Try to construct the matrix using SVDs, or Lanczos if it
                % is a sparse matrix. 
                obj = hmatrix_from_full(H, varargin{1});
            end
        end                
    end
end

