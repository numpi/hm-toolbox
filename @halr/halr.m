classdef halr
    %HALR Hierarchical matrix with arbitrary partition
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
    % H = HALR(A) constructs an HALR representation of the matrix A,
    %   which might be dense or sparse, with the default partitioning.
    %
    % To specify an arbitrary partitioning, the user has to construct a
    % quad-tree using the HALR structure itself, in the following way:
    %
    %  cluster = halr;
    %  cluster.sz [ m n ];
    %  cluster.admissible = true / false;
    %
    % which sets the size and the admissibility (i.e., if it is low-rank)
    % of a certain block. Then, it is possible to attach recursively other
    % children to it to represent lower levels, using the fields A11, A12,
    % A21, and A22. For instance, a 2-level partitioning with dense
    % diagonal blocks can be achieved as follows:
    %
    %  c = halr; c.sz = [ m n ];
    %  c.A11 = halr; c.A11.sz = [ m1 n1 ]; c.A11.admissible = false;
    %  c.A22 = halr; c.A22.sz = [ m-m1, n-n1 ]; c.A22.admissible = false;
    %  c.A21 = halr; c.A22.sz = [ m-m1, n1 ]; c.A22.admissible = true;
    %  c.A12 = halr; c.A12.sz = [ m1, n-n1 ]; c.A12.admissible = true;
    %
    % To construct an HALR with prescribed structure, one can then use
    % the command:
    %
    %  H = HALR(A, 'cluster', c);
    %
    % The same syntax can be used with all the other constructors, listed
    % here:
    %
    % H = HALR('banded', A) constructs a representation of a banded
    %     matrix given in sparse format.
    %
    % H = HALR('handle', Afun, m, n) builds an M x N matrix by evaluting
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
        function obj = halr(varargin)
            obj.admissible = false;
            
            if nargin == 0
                return;
            end
            
            % Identify the block we are working on
            progress_fcn = @(l,i,j,s) 1;
            
            % Similarly to @hodlr and @hss, we allow the user to specify a
            % custer, which needs to be an empty @halr object with just
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
            
            for j = 1 : length(varargin)
                if strcmp(varargin{j}, 'progress')                    
                    % Filter out the progress option
                    varargin = varargin([ 1 : j-1, j+1:length(varargin) ]);
                    progress_fcn = @(l,i,j,s) fprintf('Level %d / block (%d, %d): status %s\n', l, i, j, s);
                    break;
                end
            end
            
            if ischar(varargin{1})
                switch varargin{1}
                    case 'adaptive'
                        if ~isempty(H)
                            error('The adaptive constructor cannot be used with a prescribed cluster');
                        end
                        
                        if nargin <= 4                        
                            obj = halr_from_adaptive(varargin{2:4}, []);
                        else
                            obj = halr_from_adaptive(varargin{2:end});
                        end
                        
                    case 'handle'
                        obj = halr_from_aca(H, varargin{2:4}, progress_fcn);
                    case 'eye'
                        obj = halr_eye(H, varargin{2});
                    case 'banded'
                        % At the moment we just rely on the Lanczos constructor
                        obj = halr_from_full(H, varargin{2});
                    case 'low-rank'
                        obj = halr_from_low_rank(H, varargin{2}, varargin{3});
                    case 'zeros'
                        obj = halr('low-rank', zeros(varargin{2}, 0), zeros(varargin{3}, 0));
                    case 'ones'
                        obj = halr('low-rank', ones(varargin{2}, 1), ones(varargin{3}, 1));
                    otherwise
                        error('Unsupported constructor');
                end
            else
                % Try to construct the matrix using SVDs, or Lanczos if it
                % is a sparse matrix.
                obj = halr_from_full(H, varargin{1});
            end
        end
    end
end

