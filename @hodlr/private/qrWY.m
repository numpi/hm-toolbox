function [Y, T, A] = qrWY(A)
%QRWY Computes QR decomposition in compact WY representation for dense matrix.
%
% [Y,T,R] = QRWY(A) computes a QR decomposition A = Q*R of a dense matrix
%     A where Q is represented in terms of its compact WY representation:
%     Q = I - Y*T*Y', with T upper triangular and Y lower triangular.

[m,n] = size(A);
if m < n,
    error('Input matrix must have more rows than columns.');
end

% NB - hard coded block size. Adjust if needed.
nb = 10;

if n <= nb,
    Y = zeros(m,n);
    beta = zeros(n,1);
    if (m == n)
        beta(end) = 0;
        Y(m,n)=1;
        n1  = n-1; 
    else
        n1 = n;
    end
    for j = 1:n1,
        [u,beta(j)] = house(A(j:end,j));
        A(j:end,j:end) = A(j:end,j:end) - ( beta(j)*u )*( u'*A(j:end,j:end) );
        A(j+1:end,j) = 0;
        Y(j:end,j) = u;
    end
    T = zeros(n);
    for j = 1:n,
       T(1:j-1,j) = -beta(j)*T(1:j-1,1:j-1)*(Y(:,1:j-1)'*Y(:,j));
       T(j,j) = beta(j);
    end
else
    % Compute QR recursively a la Elmroth and Gustavson.
    n1 = floor(n/2); 
    j1 = n1+1;
    
    [Y1, T1, A(:, 1:n1)] = qrWY( A(:, 1:n1) );
    x = Y1*T1';
    A(:, j1:n) = A(:, j1:n) - (x*(Y1'*A(1:m, j1:n)));
   
    [Y2, T2, A(j1:m,j1:n)] = qrWY( A(j1:m,j1:n) );
    Y2 = [zeros(size(Y1,1)-size(Y2,1),size(Y2,2));Y2]; 
        
    Y = [Y1,Y2];   
    m2 = size(T1,2);
    n2= size(T2,1);
    T3 = -x'*(Y2*T2);
    
    T = [T1 T3; zeros(n2,m2) T2];
end

end

function [v,beta] = house(x)
% HOUSE
%
% Given a vector x in R^n, this computes v in R^n
% and beta in R, such that (eye - beta v*v') x = [* 0]';
%

n = length(x);
nrm = norm(x);
if nrm ~= 0
    v1 = x(1)/nrm;
    v1 = v1 + sign(v1) + ( v1==0 );
    beta = abs(v1);
    v1 = v1*nrm;
    v = [1; x(2:end) / v1];
else
    v = [1; x(2:n)];
    beta = 0;
end

end
