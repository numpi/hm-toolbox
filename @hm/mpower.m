function T = mpower(T, p)
%MPOWER Evaluate the matrix power of T.

n = size(T, 1);

if ~isscalar(p) || floor(p) ~= p
    error('Only integer powers are supported');
end

if T.sz(1) ~= T.sz(2)
    error('A^p is only possible for square matrices');
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
        T = hm('diagonal', ones(n, 1));
    otherwise
        p1 = floor(p / 2);
        p2 = p - p1;
        
        % Call this function recursively
        if p1 == p2
            T = T^p1;
            T = T*T;
        else
            T = T^p1 * T^p2;
        end
end

end

