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
