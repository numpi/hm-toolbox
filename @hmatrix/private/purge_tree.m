function H = purge_tree(H)
%PURGE_TREE Clean empty nodes in the tree of H
if isempty(H)
	% empty matrix -- no need for anything
else
	if is_leafnode(H)
    		% We are in a leafnode -- no need for anything
	else
   	 	collapsed = false;
    
    		if H.A11.sz(1) == 0 && H.A11.sz(2) == 0
        		H.A11 = [];
        		H = H.A22;
        		collapsed = true;
    		end
    
    		if ~collapsed && H.A22.sz(1) == 0 && H.A22.sz(2) == 0  
        		H.A22 = [];
        		H = H.A11;
        		collapsed = true;
    		end
    
    		if ~collapsed
        		H.A11 = purge_tree(H.A11);
        		H.A22 = purge_tree(H.A22);
                H.A12 = purge_tree(H.A12);
                H.A21 = purge_tree(H.A21);
    		else
        		H = purge_tree(H);
    		end
	end

end

