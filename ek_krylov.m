function [V, K, H, params] = ek_krylov(varargin)
%EK_KRYLOV Extended Krylov projection of a matrix A. 
%
% [V, K, H, params] = EK_KRYLOV(A, B) construct the extended Krylov 
%     subspace spanned by [B, A*B, A\B]. The matrix V is an orthogonal 
%     basis for this space, and K and H are block upper Hessenberg 
%     rectangular matrices satisfying
%
%        A * V * K = V * H                                          (1)
%
% [V, K, H, params] = EK_KRYLOV(V, K, H, PARAMS) enlarges an extended
%     Krylov subspace generated with a previous call to EK_KRYLOV by adding
%     another zero / infinity pole pair. The resulting space will satisfy
%     the same relation (1). 
%
% Note: the poles are added as [0, inf] (in this order). 

if nargin ~= 2 && nargin ~= 4
    error('Called with the wrong number of arguments');
end

if nargin == 2
    % Start to construct the extended Krylov space
    [V, K, H, params] = ek_krylov_start(varargin{:});
else
    % Enlarge the space that was previously built
    [V, K, H, params] = ek_krylov_extend(varargin{:});
end

end

function [V, K, H, params] = ek_krylov_start(A, b)

    if ~isstruct(A)
        m = size(A, 1);
        n = size(A, 2);
    
        if m ~= n
            error('The matrix A should be square');
        end

        if n ~= size(b, 1)
            error('The block vector b has wrong number of rows');
        end
    end
    
    bs = size(b, 2);
    
    % Construct a basis for the column span of b
    [V, ~] = qr(b, 0);
    
    H = zeros(bs, 0);
    K = zeros(bs, 0);
    
    [V, K, H, w] = add_zero_pole(V, K, H, A, V);
    [V, K, H, w] = add_inf_pole (V, K, H, A, w);
    
    % Save parameters for the next call
    params = struct();
    params.last = w;
    params.A = A;
end

function [V, K, H, params] = ek_krylov_extend(V, K, H, params)
    w  = params.last;
    A  = params.A;
    
    [V, K, H, w] = add_zero_pole(V, K, H, A, w);
    [V, K, H, w] = add_inf_pole (V, K, H, A, w);
    
    params.last = w;
end

%
% Utility routine that adds a zero pole to the space. The vector w is
% the continuation vector. 
%
function [V, K, H, w] = add_zero_pole(V, K, H, A, w)
    bs = size(w, 2);
    
    if isstruct(A)
        w = A.solve(1.0, 0.0, w);
    else
        w = A \ w;
    end
    
     % Perform orthogonalization with modified Gram-Schimidt
    [w, h] = mgs_orthogonalize(V, w);

    % Enlarge H and K
    H(size(H, 1) + bs, size(H, 2) + bs) = 0;
    K(size(K, 1) + bs, size(K, 2) + bs) = 0;

    K(1:end-bs, end-bs+1:end) = h;
    H(end-2*bs+1:end-bs, end-bs+1:end) = eye(bs);

    [w, r] = qr(w, 0);    
    K(end-bs+1:end, end-bs+1:end) = r;

    V = [V, w];
end

%
% Utility routine that adds an infinity pole to the space. The vector w is
% the continuation vector. 
% 
function [V, K, H, w] = add_inf_pole(V, K, H, A, w)
    bs = size(w, 2);
    
    if isstruct(A)
        w = A.multiply(1.0, 0.0, w);
    else
        w = A * w;
    end

    % Perform orthogonalization with modified Gram-Schimidt
    [w, h] = mgs_orthogonalize(V, w);

    % Enlarge H and K
    H(size(H, 1) + bs, size(H, 2) + bs) = 0;
    K(size(K, 1) + bs, size(K, 2) + bs) = 0;

    H(1:end-bs, end-bs+1:end) = h;
    K(end-2*bs+1:end-bs, end-bs+1:end) = eye(bs);

    [w, r] = qr(w, 0);    
    H(end-bs+1:end, end-bs+1:end) = r;

    V = [V, w];
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

