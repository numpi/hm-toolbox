function [U,V] = hss_offdiag(A, ul)

[Ul,Vl] = hss_generators(A.hssl);
[Ur,Vr] = hss_generators(A.hssr);

if ~exist('ul', 'var')
	ul = 'all';
end

switch ul
	case 'upper'
		U = Ul * A.Bu;
		V = Vr;
	case 'lower'
		U = Ur * A.Bl;
		V = Vl;
	case 'all'
		U = [ zeros(size(Ul,1), size(A.Bl,2)) , Ul * A.Bu ; Ur * A.Bl , zeros(size(Ur,1), size(A.Bu,2)) ];
		V = [ Vl , zeros(size(Vl,1), size(Vr,2)) ; zeros(size(Vr,1), size(Vl,2)) , Vr ];
	end
end

function [U,V] = hss_generators(A)
	if A.leafnode == 1
		U = A.U;
		V = A.V;
	else
		[Ul,Vl] = hss_generators(A.hssl);
		[Ur,Vr] = hss_generators(A.hssr);
		U = [ Ul * A.Rl ; Ur * A.Rr ];
		V = [ Vl * A.Wl ; Vr * A.Wr ];
	end
end
