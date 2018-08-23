function H = mrdivide(H1, H2)
%MRDIVIDE solve linear systems with HODLR-matrices
if isa(H2,'hm')
    H = (H2'\H1')';
elseif isscalar(H2)
    H = H1 *(1/H2);
end
