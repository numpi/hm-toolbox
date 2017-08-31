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
			B.hssl.Rl = A.hssl.Wl;
			B.hssl.Rr = A.hssl.Wr;
			B.hssl.Wl = A.hssl.Rl;
			B.hssl.Wr = A.hssl.Rr;
			B.hssr.Rl = A.hssr.Wl;
			B.hssr.Rr = A.hssr.Wr;
			B.hssr.Wl = A.hssr.Rl;
			B.hssr.Wr = A.hssr.Rr;
		else
			B.hssl.U = A.hssl.V;
			B.hssl.V = A.hssl.U;
			B.hssr.U = A.hssr.V;
			B.hssr.V = A.hssr.U;
		end
		B.hssl = hss_transpose(B.hssl);
		B.hssr = hss_transpose(B.hssr);
	end
end
