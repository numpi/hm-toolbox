function b = isreal(H)
	if is_leafnode(H)
		b = isreal(H.F);
	else
		b = isreal(H.A11) && isreal(H.A22) && isreal(H.U21) && isreal(H.V21) && isreal(H.U12) && isreal(H.V12);
	end
end
