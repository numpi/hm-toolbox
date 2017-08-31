function sz = size(H, idx)
%SIZE Size of an HSS matrix. ù

if ~isempty(H.ml)
    m = H.ml + H.mr;
else
    m = size(H.D, 1);
end

if ~isempty(H.nl)
    n = H.nl + H.nr;
else
    n = size(H.D, 2);
end

if exist('idx', 'var')
    switch idx
        case 1
            sz = m;
        case 2
            sz = n;
    end
else
    sz = [m n];
end

end

