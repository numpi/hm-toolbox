function [ U, V ] = compress_factors(Uold, Vold, nrm)
%COMPRESS_FACTORS Compress a low-rank representation U*V'.

threshold = hmatrixoption('threshold');

if isempty(Uold)
    U = Uold;
    V = Vold;
else
    [QU, RU] = qr(Uold, 0);
    [QV, RV] = qr(Vold, 0);
    
    [U,S,V] = svd(RU * RV');
    
    if ~exist('nrm', 'var') || isempty(nrm)
        switch hmatrixoption('norm')
            case 2
                nrm = S(1,1);
            case 'fro'
                nrm = norm(diag(S));
        end
    end
    
    % Estimate the rank: the code changes depending on the norm type
    t = diag(S);
    switch hmatrixoption('norm')
        case 2
            rk = sum(t > nrm * threshold);
        case 'fro'
            t = cumsum(t(end:-1:1));
            rk = sum(t > nrm * threshold);
    end
    
    U = QU * U(:,1:rk) * S(1:rk,1:rk);
    V = QV * V(:,1:rk);
end

end

