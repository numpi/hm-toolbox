function sz = size(H, idx)
%SIZE Get the size of an H-matrix.

if ~exist('idx', 'var')
    sz = H.sz;
else
    switch idx
        case 1
            sz = H.sz(1);
        case 2
            sz = H.sz(2);
        otherwise
            error('Unsupported index specified');
    end
end

end

