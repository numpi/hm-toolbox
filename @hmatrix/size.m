function varargout = size(H, idx)
%SIZE Get the size of an HMATRIX.

if nargout == 2
    varargout{1} = H.sz(1);
    varargout{2} = H.sz(2);
    return
end

if ~exist('idx', 'var')
    varargout{1} = H.sz;
else
    switch idx
        case 1
            varargout{1} = H.sz(1);
        case 2
            varargout{1} = H.sz(2);
        otherwise
            error('Unsupported index specified');
    end
end

end

