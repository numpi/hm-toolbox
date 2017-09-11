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

function B = hss_scalar_mul(a, B)
	if B.leafnode == 1		
		B.D = a * B.D;
	else
		B.Bu = a * B.Bu;
		B.Bl = a * B.Bl;
		B.hssl = hss_scalar_mul(a, B.hssl);
		B.hssr = hss_scalar_mul(a, B.hssr);
	end
end
