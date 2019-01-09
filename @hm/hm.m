classdef hm
%HM
%
% H = HM(A) constructs the HODLR representation of a matrix A by
%     approximating the off-diagonal blocks at all levels using truncated 
%     SVDs.
% 
% H = HM('low-rank', U, V) construct an HM representation of the low-rank
%     matrix U*V'. 
%
% H = HM('banded', A) constructs an HM representation of a banded matrix
%     A. The matrix A can be either sparse or dense. 
%
% H = HM('banded', A, B) can be used to specify the symmetric bandwidth B 
%     of the matrix A.
%
% H = HM('banded', A, BL, BU) specifies different lower and upper 
%     bandwidth BL and BU, respectively.
%
% H = HM('chebfun2', F, XDOM, YDOM, M, N) constructs the M x N matrix
%     containing the samplings of the bivariate function F over a uniform
%     grid of the square XDOM x YDOM. The procedure relies on separable
%     approximation of F(X,Y) as provided by the Chebfun package. 
%
% H = HM('toeplitz', C, R) constructs the Toeplitz matrix with C as first
%     column and R as first row. 
%
% H = HM('toeplitz', C, R, N) constructs the N x N Toeplitz matrix with 
%     C as first column and R as first row. if N is larger than C and 
%     R, respectively, then C and R are padded with zeros. 
    
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
                
                obj = hm_build_hm_tree(size(A, 1), size(A, 2), ...
                            hmoption('block-size'), rowcluster, ...
                            colcluster);
                
                obj = create_h_matrix(obj, varargin{1});
                return;
            end
            
            if nargin > 1
                switch varargin{1}
                    case 'low-rank'
                        obj = hm_build_hm_tree(size(varargin{2}, 1), ...
                            size(varargin{3}, 1), ...
                            hmoption('block-size'), rowcluster, ...
                            colcluster);
                        
                        obj = create_low_rank_h_matrix(obj, varargin{2:charpos - 1});
                    case 'tridiagonal'
                        obj = create_tridiagonal_h_matrix(obj, varargin{2});
                    case 'banded'
                        obj = hm_build_hm_tree(size(varargin{2}, 1), ...
                            size(varargin{2}, 2), ...
                            hmoption('block-size'), rowcluster, ...
                            colcluster);                        
                        
                        obj = create_banded_h_matrix(obj, varargin{2:charpos - 1});
                    case 'eye'
                        n = varargin{2};
                        obj = hm('diagonal', ones(n, 1), 'cluster', rowcluster);
                    case 'diagonal'
                        D = varargin{2}; D = D(:);
                        obj = hm('banded', spdiags(D, 0, ...
                            length(D), length(D)), varargin{3:end});
                    case 'chebfun2'
                        if charpos < 6
                            error('Unsufficient arguments for the chebfun2 constructor');
                        end
                        m = varargin{5};
                        if charpos == 6
                            n = m;
                        else
                            n = varargin{6};
                        end

                        obj = hm_build_hm_tree(m, n, ...
                            hmoption('block-size'), rowcluster, ...
                            colcluster); 

                        obj = create_chebfun2_h_matrix(obj, varargin{2:4}, m, n);
                    case 'toeplitz'
                        n = length(varargin{2});
                        if charpos > 4
                            n = varargin{4};
                        end
                        
                        obj = hm_build_hm_tree(n, n, ...
                            hmoption('block-size'), rowcluster, ...
                            colcluster); 
                        
                        obj = create_toeplitz_h_matrix(obj, varargin{2:charpos-1});
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
            
            if ~isempty(obj.F)
                obj.F = A;
            else
                % Get the middle point
                % mp = ceil(size(A, 1) / 2);
                m1 = obj.A11.sz(1);
                n1 = obj.A11.sz(2);
                
                obj.A11 = create_h_matrix(obj.A11, A(1:m1,1:n1));
                obj.A22 = create_h_matrix(obj.A22, A(m1+1:end,n1+1:end));
                
                [obj.U21, obj.V21] = compress_matrix(A(m1+1:end,1:n1));
                [obj.U12, obj.V12] = compress_matrix(A(1:m1,n1+1:end));
            end
        end
        
        function H = create_banded_h_matrix(obj, A, bandl, bandu)
            %CREATE_BANDED_H_MATRIX Create a banded H-matrix.
            H = obj;
            
            if ~exist('bandl', 'var')
                [bandl, bandu] = bandwidth(A);
            end
            
            if ~exist('bandu', 'var')
                bandu = bandl;
            end
            
            if ~isempty(H.F)
                H.F = full(A);
            else
                [m, n] = size(H);

                if m ~= n
                    error('The construction of banded matrices is supported only for square diagonal blocks');
                end
                
                m1 = H.A11.sz(1);
                n1 = H.A11.sz(2);
                
                if max(bandu, bandl) <= min(n - max(m1, n1), m1)
                    H.A11 = create_banded_h_matrix(H.A11, A(1:m1,1:n1), bandl, bandu);
                    H.A22 = create_banded_h_matrix(H.A22, A(m1+1:end,n1+1:end), bandl, bandu);
                    
                    H.U12 = [ zeros(m1 - bandu, bandu) ; full(A(m1 - bandu + 1:n1, m1 + 1:m1 + bandu)) ];
                    H.V12 = [ eye(bandu) ; zeros(n - n1 - bandu, bandu) ];
                    
                    H.U21 = [ full(A(m1 + 1:m1 + bandl, m1 - bandl + 1:m1)) ; zeros(m - m1 - bandl, bandl) ];
                    H.V21 = [ zeros(n1 - bandl, bandl) ; eye(bandl) ];
                else
                    H = create_h_matrix(H, full(A));
                end
            end
        end
        
        function obj = create_low_rank_h_matrix(obj, U, V)
            %CREATE_LOW_RANK_H_MATRIX Create a low rank H matrix.
            
            if ~isempty(obj.F)
                obj.F = U * V';
            else
                m1 = obj.A11.sz(1);
                n1 = obj.A11.sz(2);
                
                obj.A11 = create_low_rank_h_matrix(obj.A11, U(1:m1,:), V(1:n1,:));
                obj.A22 = create_low_rank_h_matrix(obj.A22, U(m1+1:end,:), V(n1+1:end,:));
                obj.U12 = U(1:m1,:);
                obj.V12 = V(n1+1:end,:);
                obj.U21 = U(m1+1:end,:);
                obj.V21 = V(1:n1,:);
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
        
        function obj = create_chebfun2_h_matrix(obj, fct, xdom, ydom, m, n)           
            
            % if ~exist('chebfun2') && ~exist('chebapprox2')
            %    error('Chebfun not found: did you forget to add it to the path?');
            % end
            
            x = linspace(xdom(1), xdom(2), n);
            y = linspace(ydom(1), ydom(2), m);
            
            if ~isempty(obj.F)
                obj.F = fct( ones(m,1) * x, y.' * ones(1,n) );
            else
                m1 = obj.A11.sz(1);
                n1 = obj.A11.sz(2);
                m2 = obj.A22.sz(1);
                n2 = obj.A22.sz(2);

                obj.A11 = create_chebfun2_h_matrix(obj.A11, fct, ...
                    [ x(1), x(n1) ], ...
                    [ y(1), y(m1) ], m1, n1);
                obj.A22 = create_chebfun2_h_matrix(obj.A22, fct, ...
                    [ x(n1+1), x(end) ], ...
                    [ y(m1+1), y(end) ], m2, n2);
                
                % Create the low-rank block A12 and A21
                [obj.U12, obj.V12] = chebfun2_low_rank(fct, ...
                    [ x(n1+1), x(end) ], ...
                    [ y(1), y(m1) ], ...
                    m1, n2);
                [obj.U21, obj.V21] = chebfun2_low_rank(fct, ...
                    [ x(1), x(n1) ], ...
                    [ y(m1+1), y(end) ], ...
                    m2, n1);
            end
            
            obj = compress_hmatrix(obj);
        end
        
        function obj = create_toeplitz_h_matrix(obj, am, ap, n)
            
            am = am(:).';
            ap = ap(:).';
            
            if ~exist('n', 'var')
                n = length(am);
                if n ~= length(ap)
                    error('Only square Toeplitz matrices are supported in @hm');
                end
            end
            
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
            
            if ~isempty(obj.F)
                obj.F = toeplitz(am, ap);
            else
                [mu, nu, ml, nl] = find_max_offdiag_size(obj);
                
                tol = hmoption('threshold');
                
                % Here the truncation is made relatively to the norm of the
                % off-diagonal blocks: we might want to do it relatively to
                % the norm of the big matrix.
                [tU21,S21,tV21] = lanczos_svd(@(v,trasp) toepmult_afun(...
                    [ am(nl+1:min(nl+ml, length(am))) , zeros(1, ml - min(nl+ml, length(am)) + ml) ], am(nl+1:-1:2), ...
                    ml, nl, v, trasp), ml, nl, tol);
                [tU12,S12,tV12] = lanczos_svd(@(v,trasp) toepmult_afun(...
                    ap(mu+1:-1:2), [ ap(mu+1:min(mu+nu, length(ap))), zeros(1, nu - min(mu+nu, length(ap)) + mu) ], ...
                    mu, nu, v, trasp), mu, nu, tol);
                
                tU21 = tU21 * sqrt(S21);
                tV21 = tV21 * sqrt(S21);
                tU12 = tU12 * sqrt(S12);
                tV12 = tV12 * sqrt(S12);
                
                obj = initialize_toeplitz_h_matrix(obj, am, ap, n, ...
                    tU12, tV12, tU21, tV21);
            end
            
        end
        
        function [mu, nu, ml, nl] = find_max_offdiag_size(obj)
            if isempty(obj.F)
                [mu1, nu1, ml1, nl1] = find_max_offdiag_size(obj.A11);                
                [mu2, nu2, ml2, nl2] = find_max_offdiag_size(obj.A22);
                
                mu = max([ size(obj.U12, 1), mu1, mu2 ]);
                nu = max([ size(obj.V12, 1), nu1, nu2 ]);
                ml = max([ size(obj.U21, 1), ml1, ml2 ]);
                nl = max([ size(obj.V21, 1), nl1, nl2 ]);
            else
                mu = 0;
                ml = 0;
                nu = 0;
                nl = 0;
            end
        end
        
        function obj = initialize_toeplitz_h_matrix(obj, am, ap, n, tU12, tV12, tU21, tV21)
            % obj.sz = [ n n ];
            
            if ~isempty(obj.F)
                obj.F = toeplitz(am(1:n), ap(1:n));
            else
                m1 = obj.A11.sz(1);
                n1 = obj.A11.sz(2);
                m2 = obj.A22.sz(1);
                n2 = obj.A22.sz(2);
                
                obj.U21 = tU21(1:m2, :);
                obj.V21 = tV21(end-n1+1:end,:);
                obj.U12 = tU12(end-m1+1:end,:);
                obj.V12 = tV12(1:n2, :);
                
                obj.A11 = initialize_toeplitz_h_matrix(obj.A11, am, ap, m1, tU12, tV12, tU21, tV21);
                obj.A22 = initialize_toeplitz_h_matrix(obj.A22, am, ap, n-m1, tU12, tV12, tU21, tV21);
            end
        end
        
        function obj = create_cauchy_h_matrix(obj, x, y)
            n = length(x);
            
            obj.sz = [ n n ];
            
            x = reshape(x, length(x), 1);
            y = reshape(y, length(y), 1);
            
            if n <= hmoption('block-size')
                obj.F = 1 ./ (x * ones(1, n) + ones(n, 1) * y.');
                for i = 1 : n
                    obj.F(i,i) = 0;
                end
            else
                mp = ceil(n / 2);
                np = ceil(n / 2);
                
                [obj.U21, obj.V21] = cauchy_lr(x(mp+1:end), y(1:np), hmoption('threshold'));
                [obj.U12, obj.V12] = cauchy_lr(x(1:mp), y(np+1:end), hmoption('threshold'));
                
                obj.A11 = hm();
                obj.A22 = hm();
                
                obj.A11 = create_cauchy_h_matrix(obj.A11, x(1:mp), y(1:np));
                obj.A22 = create_cauchy_h_matrix(obj.A22, x(mp+1:end), y(np+1:end));
            end
        end
        
        
    end
    
end

