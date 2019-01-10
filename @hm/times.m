function H = times(H1, H2)
%TIMES Scalar multiplication

if (isfloat(H1) && isscalar(H1)) || (isfloat(H2) && isscalar(H2)) % if one of the two is a scalar
    	H = H1 * H2;
elseif isa(H1, 'hm') && isa(H2, 'hm')
    	H = hm_hadamard_mul(H1, H2);
    	H = compress_hmatrix(H);
elseif isa(H1, 'hm') && isa(H2, 'hss')
    	H = hm_hadamard_mul(H1, hss2hm(H2));	
    	H = compress_hmatrix(H);
elseif isa(H1, 'hm') && ~isa(H2, 'hm')
	H = full(H1) .* H2;
elseif ~isa(H1, 'hm') && isa(H2, 'hm')
	H = H1 .* full(H2);
else
    error('A .* B: Unsupported operation');
end

end

