function H = solve_upper_triangular(H1, H2)
if isa(H2,'hodlr') %case of hierarchical right-hand side
    if is_leafnode(H1)
        H = H2;
        H.F = H1.F \ H2.F;
    else
        H = H2;
        %H.U21 = H1.A22\H2.U21; 
	H.U21 = solve_upper_triangular(H1.A22, H2.U21);
        %H.A22 = H1.A22\H2.A22;
	H.A22 = solve_upper_triangular(H1.A22, H2.A22);
        H.U12 = [H2.U12, H1.U12];
        H.V12 = [H2.V12, -hodlr_mtimes_dense(H.A22', H1.V12)];
        %H.U12 = H1.A11\ H.U12;
	H.U12 = solve_upper_triangular(H1.A11, H.U12);
        %H.A11 = H1.A11\hodlr_rank_update(H2.A11, -H1.U12 * (H1.V12' * (H1.A22 \ H2.U21) ) ,H2.V21);
        H.A11 = solve_upper_triangular(H1.A11, hodlr_rank_update(H2.A11, -H1.U12 * (H1.V12' * (H.U21) ) ,H2.V21));
    end
else % case of dense right-hand side
    if is_leafnode(H1)
        H = H1.F\H2;
    else
        mp = size(H1.A11,2);
        
        if ~is_leafnode(H1.A22) && isempty(H1.A22.U21)
            x2 = solve_upper_triangular(H1.A22, H2(mp+1:end, :));
        else
            x2 = H1.A22\H2(mp+1:end,:);
        end
        
        if ~is_leafnode(H1.A11) && isempty(H1.A11.U21)
            x1 = solve_upper_triangular(H1.A11, H2(1:mp,:) - H1.U12 * (H1.V12' * x2));
        else
            x1 = H1.A11\(H2(1:mp,:) - H1.U12 * (H1.V12' * x2));
        end
        
        H = [x1; x2];
    end
end
