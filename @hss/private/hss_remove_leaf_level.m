function X = hss_remove_leaf_level(X)
%


if X.leafnode == 1
    error('Cannot remove a level on a single level HSS');
else
    if X.A11.leafnode == 1 && X.A22.leafnode == 1
        was_topnode = (X.topnode == 1);
        X.topnode = 1;
        
        X.D = full(X);
        X.topnode = was_topnode;
        
        if ~was_topnode
            X.U = [ X.A11.U * X.Rl ; X.A22.U * X.Rr ];
            X.V = [ X.A11.V * X.Wl ; X.A22.V * X.Wr ];
        else
            X.Rl = []; X.Rr = []; X.Wr = []; X.Wl = [];
        end
        
        X.leafnode = 1;
        X.A11 = []; X.A22 = [];
        X.ml = []; X.nl = []; X.mr = []; X.nr = [];
    else
	if X.A11.leafnode == 0
        	X.A11 = hss_remove_leaf_level(X.A11);
	end
	if X.A22.leafnode == 0
        	X.A22 = hss_remove_leaf_level(X.A22);
	end
    end
end


end

