function H = times(H1, H2)
%TIMES Scalar multiplication

if (isfloat(H1) && isscalar(H1)) || (isfloat(H2) && isscalar(H2)) % if one of the two is a scalar
    H = H1 * H2;
elseif isa(H1, 'hodlr') && isa(H2, 'hodlr')
    H = hodlr_hadamard_mul(H1, H2);
    H = compress_hodlr(H);
elseif isa(H1, 'hodlr') && isa(H2, 'hss')
    H = hodlr_hadamard_mul(H1, hss2hodlr(H2));
    H = compress_hodlr(H);
elseif isa(H1, 'hodlr') && ~isa(H2, 'hodlr')
    H = full(H1) .* H2;
elseif ~isa(H1, 'hodlr') && isa(H2, 'hodlr')
    H = H1 .* full(H2);
else
    error('A .* B: Unsupported operation');
end

end

