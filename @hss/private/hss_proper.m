function B = hss_proper(A)
	% Forward stage
	if (A.leafnode==1)
		B = A;
		return; % It means that we have called the compression on a HSS matrix with only the root, so nothing to do
	end
	if (A.hssl.leafnode==1)
		[Ut, Pt] = qr(A.hssr.U,0);
		[Vt, Qt] = qr(A.hssr.V,0);
		A.hssr.U = Ut;
		A.hssr.V = Vt;
		if (A.topnode==0)
			A.Rr = Pt * A.Rr;
			A.Wr = Qt * A.Wr;
		end
		A.Bl = Pt * A.Bl;  
		A.Bu = A.Bu * Qt';
		[Ut, Pt] = qr(A.hssl.U,0);
		[Vt, Qt] = qr(A.hssl.V,0);
		A.hssl.U = Ut;
		A.hssl.V = Vt;
		if (A.topnode==0)
			A.Rl = Pt * A.Rl;
			A.Wl = Qt * A.Wl;
		end
		A.Bl = A.Bl * Qt';  
		A.Bu = Pt * A.Bu;
		B = A;
	else
		A.hssl = hss_proper(A.hssl);
		A.hssr = hss_proper(A.hssr);

		[Ut, Pt] = qr([A.hssl.Rl;A.hssl.Rr],0);
		[Vt, Qt] = qr([A.hssl.Wl;A.hssl.Wr],0);
		j = size(A.hssl.Rl,1);
		A.hssl.Rl = Ut(1:j,:);
		A.hssl.Rr = Ut(j+1:end,:);
		j = size(A.hssl.Wl,1);
		A.hssl.Wl = Vt(1:j,:);
		A.hssl.Wr = Vt(j+1:end,:);
		if (A.topnode==0)
			A.Rl = Pt * A.Rl;
			A.Wl = Qt * A.Wl;
		end
		A.Bl = A.Bl * Qt';  
		A.Bu = Pt * A.Bu;
		
		[Ut, Pt] = qr([A.hssr.Rl;A.hssr.Rr],0);
		[Vt, Qt] = qr([A.hssr.Wl;A.hssr.Wr],0);
		j = size(A.hssr.Rl,1);
		A.hssr.Rl = Ut(1:j,:);
		A.hssr.Rr = Ut(j+1:end,:);
		j = size(A.hssr.Wl,1);
		A.hssr.Wl = Vt(1:j,:);
		A.hssr.Wr = Vt(j+1:end,:);
		if (A.topnode==0)
			A.Rr = Pt * A.Rr;
			A.Wr = Qt * A.Wr;
		end
		A.Bl = Pt * A.Bl;  
		A.Bu = A.Bu * Qt';
		B = A;
	end
end

