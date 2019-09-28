function [V, K, H, params] = rk_krylov(varargin)
%RK_KRYLOV Rational (block) Krylov projection of a matrix A. 
%
% [V, K, H, PARAMS] = RK_KRYLOV(A, U, POLES) construct the rational Krylov 
%     subspace spanned by 
% 
%      [U, (A - POLES(1) * I) \ U, ..., (A - POLES(L)*I) \ U], 
%
%     where L is the length of the vector V. The matrix V is an orthogonal 
%     basis for this space, and K and H are block upper Hessenberg 
%     rectangular matrices satisfying
%
%        A * V * K = V * H                                          (1)
%
% [V, K, H, PARAMS] = RK_KRYLOV(A, B, U, POLES) constructs the analogous
%     space spanned by the pencil (A, B), i.e., 
%
%      [U, (A - POLES(1) * B) \ U, ..., (A - POLES(L) * B) \ U], 
%
%     In this case, the matrix relation is A*V*K = B*V*H (2). 
% 
% [V, K, H, PARAMS] = RK_KRYLOV(A, B, V, K, H, POLES, PARAMS) enlarges a 
%     rational Krylov subspace generated with a previous call to RK_KRYLOV 
%     by adding the new poles in the vector POLES. The resulting space will
%     satisfy the same relation (2). 
%
% The parameters are inspired by RKTOOLBOX, with the idea that the latter 
% (which is much more complete) may be used as a drop-in replacement for 
% this file. 

if nargin ~= 3 && nargin ~= 4 && nargin ~= 6 && nargin ~= 7
    error('Unsupported number of parameters');
end

A = varargin{1};

if ~isstruct(A) && ( nargin == 3 || nargin == 6 )
    B = eye(size(A), 'like', A);
else
    B = varargin{2};
end

if ~isstruct(A)
    A = rk_struct(A, B);
end

if nargin == 3 || nargin == 4
    U = varargin{nargin - 1};
    poles = varargin{nargin};
    bs = size(U, 2);
    
    [V, K, H] = rk_krylov_start(U);
    
    params = struct('bs', bs);
else
    j = 1;
    if nargin == 7
        j = 2;
    end
    
    V = varargin{j+1};
    K = varargin{j+2};
    H = varargin{j+3};
    poles = varargin{j+4};
    params = varargin{j+5};
    bs = params.bs;
end

for j = 1 : length(poles)
    [V, K, H] = rk_add_pole(A, V, K, H, poles(j), bs);
end

end

function [V, K, H] = rk_krylov_start(U)
    bs = size(U, 2);
    K = zeros(bs, 0);
    H = zeros(bs, 0);
    [V, ~] = qr(U, 0);
end

function [V, K, H] = rk_add_pole(A, V, K, H, pole, bs)
    % FIXME: The choice of the continuation vector is not properly
    % implemented yet -- this part should be improved. 
    c = kron(continuation_vector(K, H, pole, bs), eye(bs));
    
    % Construct the new elements in the bottom-right corners of K and H
    if abs(pole) > 1
        kh = 1 / pole;
        hh = 1;

        % Solve the shifted linear system
        % w = ( kh * A - hh * eye(n, 'like', A) ) \ (A * V * c);
        w = correction( A.solve(kh, hh, A.multiply(1, 0, cqt(V * c))) );
    else
        kh = 1;
        hh = pole;
        
        % Solve the shifted linear system
        % w = ( kh * A - hh * eye(n, 'like', A) ) \ (V * c);
        w = correction( A.solve(kh, hh, -A.multiply(0, 1, cqt(V * c))) );
    end
        
    % We might need to add more rows to the basis
    if size(w, 1) > size(V, 1)
        V(size(w, 1), 1) = 0;
    end
    if size(V, 1) > size(w, 1)
        w(size(V, 1), 1) = 0;
    end
    
    % Perform re-orthogonalization
    [w, b] = mgs_orthogonalize(V, w);
    
    % Normalize w
    [w, s] = qr(w, 0);
    
    % Extend the matrices
    if abs(pole) > 1
        K = [ K , c - kh * b ; zeros(bs, size(K, 2)) , -s * kh ];
        H = [ H , -hh * b ; zeros(bs, size(H, 2)) , -s * hh ];
    else
        K = [ K , -kh * b ; zeros(bs, size(K, 2)) , -s * kh ];
        H = [ H , -c - hh * b ; zeros(bs, size(H, 2)) , -s * hh ];
    end
    
    V = [ V, w ];
end

function c = continuation_vector(K, H, pole, bs)
    % Find the list of poles we had

    c = ones(size(H, 1) / bs, 1);
end

function p = poles_pencil(K, H, bs)
    j = 0;
    p = zeros(1, size(K, 2) / bs);
    
    % while j < size(K, 2)
    %    l = eig( 
    %    p(j / bs + 1) = 
    %end
end

%
% Modified Gram-Schmidt orthogonalization procedure. 
%
% Suggested improvements: work with block-size matrix vector products to
% get BLAS3 speeds. 
% 
function [w, h] = mgs_orthogonalize(V, w)
    h = zeros(size(V, 2), size(w, 2));
    for j = 1 : size(V, 2)
        h(j,:) = (V(:,j)' * w);
        w = w - V(:,j) * (V(:,j)' * w);
    end
end