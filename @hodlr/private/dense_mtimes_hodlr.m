function w = dense_mtimes_hodlr(w, H)
%DENSE_MTIMES_HODLRATRIX Dense product w * H

if is_leafnode(H)
    w = w * H.F;
else
    mp = H.A11.sz(1);
    
    w = [ dense_mtimes_hodlr(w(:,1:mp), H.A11) + ...
        (w(:,mp+1:end) * H.U21) * H.V21.' , ...
        (w(:,1:mp) * H.U12) * H.V12.' + ...
        dense_mtimes_hodlr(w(:,mp+1:end), H.A22) ];
end


end

