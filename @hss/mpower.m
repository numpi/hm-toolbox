function T = mpower(T, p)
%MPOWER Evaluate the matrix power of T.

n = size(T, 1);

if ~isscalar(p) || floor(p) ~= p
    error('Only integer powers are supported');
end

if ~check_cluster_equality(T)
    error('A^p is only supported for square matrices with square diagonal blocks');
end

if p < 0
    T = inv(T)^(-p);
    return;
end

switch p
    case 2
        T = T * T;
    case 1
        % Nothing to do
    case 0
        [r,c] = cluster(T);
        T = hss('eye', n, ones(n, 1), 'cluster', r);
    otherwise
        b = dec2bin(p);
        if b(end) == '1',
            P = T;
        else
            [r,c] = cluster(T);
            P = hss('eye', n, ones(n, 1), 'cluster', r);
        end
        for j = length(b)-1:-1:1,
            T = T*T;
            if b(j) == '1', P = P*T; end
        end
        T = P;
end

end

