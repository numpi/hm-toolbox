function C = mtimes(A,B)
if isa(B,'hss')
	if isa(A,'hss')
		error('unsupported')
	else
		C = hss_mul(B.',A.').';
	end	
else
	C = hss_mul(A,B);
end
end
