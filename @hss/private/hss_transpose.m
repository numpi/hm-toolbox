function B = hss_transpose(A)
	B = A;
	if B.leafnode == 1;
		B.D = B.D.';
	else
		B.ml = A.nl;
		B.nl = A.ml;
		B.mr = A.nr;
		B.nr = A.mr;
		B.Bu = A.Bl.';
		B.Bl = A.Bu.';
		if B.hssl.leafnode == 0
			B.hssl.Rl = conj(A.hssl.Wl);
			B.hssl.Rr = conj(A.hssl.Wr);
			B.hssl.Wl = conj(A.hssl.Rl);
			B.hssl.Wr = conj(A.hssl.Rr);
			B.hssr.Rl = conj(A.hssr.Wl);
			B.hssr.Rr = conj(A.hssr.Wr);
			B.hssr.Wl = conj(A.hssr.Rl);
			B.hssr.Wr = conj(A.hssr.Rr);
		else
			B.hssl.U = conj(A.hssl.V);
			B.hssl.V = conj(A.hssl.U);
			B.hssr.U = conj(A.hssr.V);
			B.hssr.V = conj(A.hssr.U);
		end
		B.hssl = hss_transpose(B.hssl);
		B.hssr = hss_transpose(B.hssr);
	end
end
