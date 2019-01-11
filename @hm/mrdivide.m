function H = mrdivide(H1, H2)
%MRDIVIDE solve linear systems with HODLR-matrices
if isa(H2, 'hm')
    if size(H2, 1) ~= size(H2, 2)
        error('A / B: Matrix B is not square');
    end
    
    if size(H1, 2) ~= size(H2, 1)
        error('A / B: Dimension mismatch');
    end
    H = (H2'\H1')';
elseif isscalar(H2)
    H = H1 *(1/H2);
elseif isa(H2, 'hss')
    H = H1/hss2hm(H2);
else
    H = full(H1)/H2;
end
