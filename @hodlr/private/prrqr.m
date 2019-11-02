function [Q, R, p] = prrqr(A, tol, debug)
%PRRQR Rank-Revealing QR factorization with early termination.
%
% [Q, R, P] = RRQR(A, TOL) computes a rank-revealing QR factorization of
%     the matrix A, with early termination in case A is low-rank up to some
%     (relative) threshold TOL.
%
%     The stopping criterion is checked against the norms of the colums.
%
% [Q, R, P] = RRQR(A, TOL, DEBUG) specifies the debug flag. If DEBUG is set
%     to TRUE, then the method stops details about the columns norms as it
%     proceeds. The default is DEBUG = false.

if ~exist('debug', 'var')
    debug = false;
end

[m, n] = size(A);

% Compute the vector of norms
v = zeros(1, n);
for j = 1 : n
    v(j) = norm(A(:,j))^2;
end

% We will use this cell array to store the Householder reflectors
Pu = {}; Pb = {};

p = 1 : n;

for j = 1 : min(m, n - 1)
    % Find the maximum norm vector, and move it to the fron
    [mx, l] = max(v(j:end));
    v(l + j - 1) = v(j);
    v(j) = mx;
    
    if debug
        fprintf('Column norm: %e, Max col norm: %e, Threshold: %e\n', ...
            sqrt(mx), sqrt(v(1)), sqrt(v(1)) * tol * sqrt(n));
    end
    
    if sqrt(mx) < tol * sqrt(v(1)) / sqrt(j)
        R = A(1:j-1, :);
        
        Q = eye(m, j-1);
        for jj = j-1 : -1 : 1
            Q(jj:end,:) = Q(jj:end,:) - Pb{jj} * Pu{jj} * (Pu{jj}' * Q(jj:end,:));
        end
        
        if nargout == 2
            ip = zeros(1, n);
            ip(p) = 1:n;
            R = R(:,ip);
        end
        
        return;
    end
    
    % Memorize the permutation of columns
    t = p(j+l-1);
    p(j+l-1) = p(j);
    p(j) = t;
    
    % Update A as well
    w = A(:,l+j-1); A(:,l+j-1) = A(:,j);
    A(:,j) = w; w = w(j:end);
    
    % Compute the Householder reflector that maps the vector to a multiple
    % of ej
    [b, u] = householder_reflector(w);
    
    % Apply it to A
    A(j:end, j:end) = A(j:end,j:end) - b * u * (u' * A(j:end,j:end));
    
    % ...and to Q -- but this has to be delayed because we do not (yet)
    % know the number of columns that are necessary in Q. So we store the
    % Householder vectors and scalars instead.
    Pu = [ Pu, u ];
    Pb = [ Pb, b ];
    
    % Downdate the norms -- we selectively recompute them to make sure
    % there is not an excess of cancellation in this operation.
    for k = j + 1 : n
        if v(k) < A(j,k) * 1.05
            v(k) = norm(A(j+1:end,k))^2;
        else
            v(k) = v(k) - A(j,k)^2;
        end
    end
end

R = A;

Q = eye(n);
for j = n - 1 : -1 : 1
    Q(j:end,:) = Q(j:end,:) - Pb{j} * Pu{j} * (Pu{j}' * Q(j:end,:));
end

if nargout == 2
    ip = zeros(1, n);
    ip(p) = 1:n;
    R = R(:,ip);
end

end

function [b, u] = householder_reflector(w)
u = w;
b = norm(u);
u(1) = u(1) + b * sign(u(1));
b = 2 / dot(u, u);
end
