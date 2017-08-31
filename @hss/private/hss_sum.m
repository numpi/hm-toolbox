function C = hss_sum(A,B)
tol = hssoption('threshold');
C = hss_sum_ric(A,B);
C = hss_compress(C,tol);
end

function C = hss_sum_ric(A,B)
	if  A.topnode ~= B.topnode || A.leafnode ~= B.leafnode
		error('HSS_SUM:: the two hss matrices have not compatible partitioning')
	else
		C.topnode = A.topnode;
		C.leafnode = A.leafnode;
		if C.leafnode == 0
			if A.ml ~= B.ml || A.nl ~= B.nl ||A.mr ~= B.mr ||A.nr ~= B.nr 
				error('HSS_SUM:: the two hss matrices have not compatible partitioning')
			end
			C.ml = A.ml;
			C.nl = A.nl;
			C.mr = A.mr;
			C.nr = A.nr;
		end
	end

	if C.leafnode == 1
		C.D = A.D + B.D;
		C.U = [A.U, B.U];
		C.V = [A.V, B.V];
	else 
		C.Bu = blkdiag(A.Bu, B.Bu);
		C.Bl = blkdiag(A.Bl, B.Bl);
		if C.topnode == 0
			C.Rl = blkdiag(A.Rl, B.Rl);
			C.Rr = blkdiag(A.Rr, B.Rr);		
			C.Wl = blkdiag(A.Wl, B.Wl);
			C.Wr = blkdiag(A.Wr, B.Wr);
		end
		C.hssl = hss_sum_ric(A.hssl, B.hssl);
		C.hssr = hss_sum_ric(A.hssr, B.hssr);	
	end
end
