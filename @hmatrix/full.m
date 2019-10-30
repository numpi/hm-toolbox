function M = full(H)
%FULL Construct the dense version of H.

if is_leafnode(H)
	if H.admissible
    		M = H.U * H.V';
	else
		M = H.F;
	end
else
    M = [ full(H.A11) , full(H.A12) ; ...
        full(H.A21) , full(H.A22) ];
end


end

