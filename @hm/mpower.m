function H = mpower(H, p)
%MPOWER Evaluate the matrix power of H.

n = size(H, 1);

if ~isscalar(p) || floor(p) ~= p
    error('Only integer powers are supported');
end

if ~check_cluster_equality(H)
    error('A^p is only supported for square matrices with square diagonal blocks');
end

if p < 0
    H = inv(H)^(-p);
    return;
end

switch p
    case 2
        H = H * H;
    case 1
        % Nothing to do
    case 0
        [r,c] = cluster(H);
        H = hm('eye', n, ones(n, 1), 'cluster', r);
    otherwise
        b = dec2bin(p);
        if b(end) == '1',
            P = H;
        else
            [r,c] = cluster(H);
            P = hm('eye', n, ones(n, 1), 'cluster', r);
        end
        for j = length(b)-1:-1:1,
            H = H * H;
            if b(j) == '1', P = P * H; end
        end
        H = P;
end

end

