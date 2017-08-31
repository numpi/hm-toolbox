function [D,U,R,B,W,V,tr] = hss2xia(A)
	if A.leafnode == 1
		D = { A.D };
		U = { A.U };
		V = { A.V };
		W = {};
		R = {};
		tr = [ 0 ];
		B = {};
	else
		[Dl,Ul,Rl,Bl,Wl,Vl,trl] = hss2xia(A.hssl);
		[Dr,Ur,Rr,Br,Wr,Vr,trr] = hss2xia(A.hssr);

		[D,U,V,tr] = mergetree(Dl,Ul,Vl,trl,Dr,Ur,Vr,trr);

		B = { Bl{:}, A.Bu, Br{:}, A.Bl };

		if A.topnode == 1
			W = { Wl{:}, [], Wr{:} };
			R = { Rl{:}, [], Rr{:} };
		else
			W = { Wl{:}, A.Wl, Wr{:}, A.Wr };
			R = { Rl{:}, A.Rl, Rr{:}, A.Rr };
		end
	end	
end

function [D,U,V,tr] = mergetree(Dl,Ul,Vl,trl,Dr,Ur,Vr,trr)
	D = { Dl{:}, Dr{:}, [] };
	U = { Ul{:}, Ur{:}, [] };
	V = { Vl{:}, Vr{:}, [] };

	s = length(trl) + length(trr) + 1;

	trl(end) = s;
	trr(end) = s - length(trl);
	tr = [ trl, trr + length(trl), 0 ];
end
