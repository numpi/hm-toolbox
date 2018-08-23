function M = full(H)
%FULL Construct the dense version of H.

if ~isempty(H.F)
    M = H.F;
else
    M = [ full(H.A11) , H.U12 * H.V12' ; ...
        H.U21 * H.V21' , full(H.A22) ];
end


end

