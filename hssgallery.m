function varargout = hssgallery(name, varargin)
%HSSGALLERY Generate examples of H-matrices.
%
% H = HSSGALLERY('laplacian', n) generates a n x n discretization of the
% Laplacian operator on [-1, 1] using a grid of n+1 points, assuming
% Dirichlet boundary conditions.
%
% H = HSSGALLERY('haber', n) generates a 6n x 6n discretization for the
% Heat diffusion problem described in [1]. 
%
% H = HSSGALLERY('rand', n, k) generates an n x n random HODLR matrix
% with HODLR rank equal to K. 
%
% [1] Haber, Aleksandar, and Michel Verhaegen. "Sparse solution of the 
%     Lyapunov equation for large-scale interconnected systems." 
%     Automatica 73 (2016): 256-268.

switch name
    case 'rand'
        if isempty(varargin)
            n = 2048;
        else
            n = varargin{1};
        end
        
        if length(varargin) < 2
            k = 1;
        else
            k = varargin{2};
        end
        
        varargout{1} = create_random_hss(n, k);
        varargout{1}.topnode = 1;
        
    case 'laplacian'
        if isempty(varargin)
            n = 2048;
        else
            n = varargin{1};
        end
        
        varargout{1} = hss('banded', ...
            spdiags(ones(n,1) * [ -1 2 -1 ], -1:1, n, n), 1, 1);
        
    case 'haber'
        if isempty('varargin')
            n = 2048;
        else
            n = varargin{1};
        end
        
        m = 6 * n;
        
        % We first generate the matrices A and P as sparse matrices.
        A = sparse(m, m);
        P = sparse(m, m);
        
        a = -1.36;
        e = 0.34;
        
        % We build the diagonal blocks
        for i = 1 : n
            range = ((i-1)*6+1) : i*6;
            A(range, range) = spdiags(ones(6,1) * [ e, a, e ], -1:1, 6, 6);
            P(range, range) = -.2 * ones(6) - .8 * eye(6);
        end
        
        % And now super and sub-diagonals
        for i = 1 : n - 1
            range = ((i-1)*6+1) : i*6;
            
            A(range, range + 6) = e * speye(6);
            A(range + 6, range) = e * speye(6);
            
            P(range, range + 6) = -.1 * ones(6);
            P(range + 6, range) = -.1 * ones(6);
        end
        
        varargout{1} = hss('banded', A, 6, 6);
        varargout{2} = compress(hss('banded', P, 11, 11));
        varargout{3} = A;
        varargout{4} = P;
    otherwise
        error('Unsupported example name');
end

    function H = create_random_hss(n, k)
        H = hss();
        
        if n <= hssoption('block-size')
            H.D = rand(n);
            [H.U, ~] = qr(rand(n, k), 0);
            [H.V, ~] = qr(rand(n, k), 0);
            
            H.leafnode = 1;
        else
            H.leafnode = 0;
            H.B12 = randn(k);
            H.B21 = randn(k);
            
            H.Rr = randn(k);
            H.Rl = randn(k);
            H.Wl = randn(k);
            H.Wr = randn(k);
            
            H.ml = ceil(n / 2);
            H.nl = H.ml;
            H.mr = n - H.ml;
            H.nr = n - H.nl;
            
            H.A11 = create_random_hss(H.ml, k);
            H.A22 = create_random_hss(H.mr, k);
        end
        
        H.topnode  = 0;
    end


end

