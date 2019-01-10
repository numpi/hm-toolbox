function H = times(H1, H2)

if (isfloat(H1) && isscalar(H1)) || (isfloat(H2) && isscalar(H2))
    H = H1 * H2;
elseif isa(H1, 'hss') && isa(H2, 'hss')
    	H = hss_hadamard_mul(H1, H2);
elseif isa(H1, 'hss') && isa(H2, 'hm')
    	H = hss2hm(H1) .* H2;	
elseif isa(H1, 'hss') && ~isa(H2, 'hm')
	H = full(H1) .* H2;
elseif ~isa(H1, 'hss') && isa(H2, 'hss')
	H = H1 .* full(H2);
else
    error('A .* B: Unsupported case');
end
end
