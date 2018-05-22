function [U, V] = cauchy_lr(x, y, tol)
%CAUCHY Construct a low-rank representation of C= 1 / (x(i) + y(j))
%
% [U, V] = CAUCHY(X) constructs a low-rank representation of the matrix
%     with entries 1 / (x(i) + y(j)).
%
% [U, V] = CAUCHY(X, TOL) truncates with the given tolerance. 

if isempty(x)
    U = [];
    V = zeros(length(y), 0);
    return;
end

if isempty(y)
    U = zeros(length(x), 0);
    V = [];
    return;
end

% Compression enabled?
compress = true;

% Minimum block-size
bs = 128;

% Tolerance
if ~exist('tol', 'var')
    tol = 1e-8;
end

nx = length(x);
ny = length(y);

x = reshape(x, nx, 1);
y = reshape(y, ny, 1);

% fprintf(' => Cauchy(%d, %d)\n', nx, ny);

if max(nx, ny) <= bs
    C = 1 ./ (x + y.');
    
    [U, S, V] = svd(C, 'econ');
    
    % Estimate the numerical rank
    rk = sum(diag(S) > S(1,1) * tol);
    
    U = U(:, 1:rk) * sqrt(S(1:rk,1:rk));
    V = V(:, 1:rk) * sqrt(S(1:rk,1:rk));
else
    xm = mean(x);
    ym = mean(y);
    
    if max(x - xm) < .5 * min(xm + y)
        % Perform a Taylor expansion in the variable x around the point xm        
        k = -floor(log(tol) / log(min(xm+y) / max(abs(x-xm))));
        scl = 1 ./ max(x - xm);
        
        U = ones(nx, k + 1);
        V = ones(ny,  k + 1);
        
        V(:,1) = 1 ./ (y + xm);
        
        for j = 1 : k
            U(:,j+1) = -U(:,j) .* (x - xm) .* scl;
            V(:,j+1) = V(:,j) ./ (y + xm) ./ scl;
        end
    elseif max(y - ym) < .5 * min(ym + x)
        % Perform a Taylor expansion in the variable y around the point ym
        k = -floor(log(tol) / log(min(ym+x) / max(abs(y-ym))));
        scl = 1 ./ max(y - ym);
        
        U = ones(nx, k + 1);
        V = ones(ny, k + 1);
        
        U(:,1) = 1 ./ (x + ym);
        
        for j = 1 : k
            U(:,j+1) = -U(:,j) ./ (x + ym) .* scl;
            V(:,j+1) = V(:,j) .* (y - ym) ./ scl;
        end
    else
        % If no expansion is possible with a fast decay, further subdivided
        % the domain, until that condition is met. 
        if nx > bs
            mpx = floor(nx / 2);
        else
            mpx = nx;
        end
        
        if ny > bs
            mpy = floor(ny / 2);
        else
            mpy = ny;
        end

        [U11, V11] = cauchy_lr(x(1:mpx), y(1:mpy));
        [U, V2]    = cauchy_lr(x, y(mpy+1:end));
        [U21, V21] = cauchy_lr(x(mpx+1:end), y(1:mpy));

        % Reconstruct the global low-rank representation
        U = [ ...
              [ U11 ; zeros(nx-mpx, size(U11,2)) ], ...
              U, ...
              [ zeros(mpx, size(U21, 2)) ; U21 ] ...
             ];

         V = [ ...
               [ V11 ; zeros(ny-mpy, size(V11,2)) ], ...
               [ zeros(mpy,size(V2,2)) ; V2 ], ...
               [ V21 ; zeros(ny-mpy,size(V21,2)) ] ...
             ];

         % Recompression
         if compress
             [QU, RU] = qr(U, 0);
             [QV, RV] = qr(V, 0);

             [U, S, V] = svd(RU * RV', 'econ');

             rk = sum(diag(S) > tol * S(1,1));

             U = QU * U(:,1:rk) * sqrt(S(1:rk,1:rk));
             V = QV * V(:,1:rk) * sqrt(S(1:rk,1:rk));
         end
    end
end
         
end