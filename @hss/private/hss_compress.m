% function that performs the compression of an HSS matrix
function B = hss_compress(A,tol)
	B = forward_stage(A);
	B = backward_stage(B,tol);
end

function B = forward_stage(A)
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
		A.hssl = forward_stage(A.hssl);
		A.hssr = forward_stage(A.hssr);

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

function A = backward_stage(A,tol,S,T)
	if(A.leafnode==1)
		return
	end
	if (A.topnode == 1)
		[U,Su,V] = tsvd(A.Bu, tol);
		A.Bu = Su; Tl = Su.';
		if A.hssl.leafnode ==0
			A.hssl.Rl = A.hssl.Rl * U;
			A.hssl.Rr = A.hssl.Rr * U;
			A.hssr.Wl = A.hssr.Wl * V;
			A.hssr.Wr = A.hssr.Wr * V;
		else
			A.hssl.U = A.hssl.U * U;
			A.hssr.V = A.hssr.V * V;
		end
		[U,Sl,V] = tsvd(A.Bl, tol);
		A.Bl = Sl; Tu = Sl.';
		if A.hssl.leafnode ==0
			A.hssr.Rl = A.hssr.Rl * U;
			A.hssr.Rr = A.hssr.Rr * U;
			A.hssl.Wl = A.hssl.Wl * V;
			A.hssl.Wr = A.hssl.Wr * V;
			A.hssl = backward_stage(A.hssl, tol, Su, Tu);
			A.hssr = backward_stage(A.hssr, tol, Sl, Tl);
		else
			A.hssr.U = A.hssr.U * U;
			A.hssl.V = A.hssl.V * V;
		end
	else
		Su = [A.Bu, A.Rl * S];		
		Tl = [A.Bu.', A.Wr * T]; % possibile inghippo
		[Us,Su,Vs] = tsvd(Su, tol);
		[Ut,Tl,Vt] = tsvd(Tl, tol);
		k = size(A.Bu,2);
		A.Bu = Su * Vs(1:k, :)' * Ut;
		A.Rl = Us' * A.Rl; 
		A.Wr = Ut' * A.Wr;
		if A.hssl.leafnode == 0
			A.hssl.Rl = A.hssl.Rl * Us;
			A.hssl.Rr = A.hssl.Rr * Us;
			A.hssr.Wl = A.hssr.Wl * Ut;
			A.hssr.Wr = A.hssr.Wr * Ut;
		else
			A.hssl.U = A.hssl.U * Us;
			A.hssr.V = A.hssr.V * Ut;
		end

		Sl = [A.Bl, A.Rr * S];		
		Tu = [A.Bl.', A.Wl * T]; % possibile inghippo
		[Us,Sl,Vs] = tsvd(Sl, tol);
		[Ut,Tu,Vt] = tsvd(Tu, tol);
		k = size(A.Bl,2);
		A.Bl = Sl * Vs(1:k, :)' * Ut;
		A.Rr = Us' * A.Rr; 
		A.Wl = Ut' * A.Wl;
		if A.hssr.leafnode == 0
			A.hssr.Rl = A.hssr.Rl * Us;
			A.hssr.Rr = A.hssr.Rr * Us;
			A.hssl.Wl = A.hssl.Wl * Ut;
			A.hssl.Wr = A.hssl.Wr * Ut;
			A.hssl = backward_stage(A.hssl, tol, Su, Tu);
			A.hssr = backward_stage(A.hssr, tol, Sl, Tl);
		else		

			A.hssr.U = A.hssr.U * Us;
			A.hssl.V = A.hssl.V * Ut;
		end
	end
end
