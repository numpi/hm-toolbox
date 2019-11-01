function [ U, V ] = compress_factors(Uold, Vold, nrm)
%COMPRESS_FACTORS Compress a low-rank representation U*V'.

threshold = hssoption('threshold');

if isempty(Uold)
    U = Uold;
    V = Vold;
else
    [QU, RU] = qr(Uold, 0);
    [QV, RV] = qr(Vold, 0);
    
    [U,S,V] = svd(RU * RV');
    
    if ~exist('nrm', 'var')
        nrm = S(1,1);
    end
    
    %rk = sum(diag(S) > threshold);
    %rk = min(sum(diag(S) > S(1,1) * threshold),50);
    
    % 2 norm of the residuals after truncation
    t = diag(S);
    % t = cumsum(t(end:-1:1));
    
    % Numerical rank of the outer product
    rk = sum(t > nrm * threshold);
    
    U = QU * U(:,1:rk) * S(1:rk,1:rk);
    V = QV * V(:,1:rk);
end

end

