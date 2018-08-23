function X = lyap(varargin)
%LYAP Solve a Lyapunov or Sylvester equation.
%
% X = LYAP(A, C) solves the Lyapunov equation AX + XA' + C = 0.
%
% X = LYAP(A, C, SA), where SA is a sparse version of A, solves the same
%     equation but exploiting the sparsity of A to build the rational
%     Krylov subspace used for the approximation.
%
% X = LYAP(A, B, C) solves the Sylvester equation AX + XB' + C = 0.
%
% X = LYAP(A, B, C, SA, SB) solves the same equation but again relying on
%     the sparsity of A and B, represented in SA and SB.

if nargin == 2
    % HSS Lyapunov case
    X = hss_dac_lyap(varargin{1}, varargin{1}', varargin{2});
elseif nargin == 3
    % Either HSS Sylvester or Sparse Lyapunov.
    if issparse(varargin{3})
        nrmA = 1 / norm(varargin{1});
        
        A  = varargin{1} * nrmA;
        C  = varargin{2} * nrmA;
        sA = varargin{3} * nrmA;
        
        X = hss_sparse_dac_lyap(A, C, sA);
    else
        X = hss_dac_lyap(varargin{1}, varargin{2}, varargin{3});
    end
elseif nargin == 5
    % Sparse Sylvester equation
    nrm = min(1 / norm(varargin{1}), 1 / norm(varargin{2}));
    
    A = varargin{1} * nrm;
    B = varargin{2} * nrm;
    C = varargin{3} * nrm;
    
    sA = varargin{4} * nrm;
    sB = varargin{5} * nrm;
    
    X = hss_sparse_dac_sylv(A, B, C, sA, sB);
else
    error('Wrong number of parameters passed to the function lyap');
end

end
