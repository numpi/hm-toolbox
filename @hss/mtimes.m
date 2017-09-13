function C = mtimes(A,B)
if isa(B,'hss')
	if isa(A,'hss')
		C = HSS_HSS_Multiply(A, B);
	elseif isscalar(A)
		C = hss_scalar_mul(A,B);
	else
		C = hss_mul(B.',A.').';
	end	
else
	if isscalar(B)
		C = hss_scalar_mul(B,A);
	else
		C = hss_mul(A,B);
	end
end
end


