function H = mpower(H, p)
%MPOWER Evaluate the matrix power of H.

n = size(H, 1);

if ~isscalar(p) || floor(p) ~= p
    error('Only integer powers are supported');
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
        H = halr('eye', n);
    otherwise
        b = dec2bin(p);
        if b(end) == '1'
            P = H;
        else
            P = halr('eye', n);
        end
        for j = length(b)-1:-1:1
            H = H * H;
            if b(j) == '1', P = P * H; end
        end
        H = P;
end

end

