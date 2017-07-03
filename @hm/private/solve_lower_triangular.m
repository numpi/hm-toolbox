function H = solve_lower_triangular(H1,H2)
if isa(H2,'hm') %case of hierarchical right-hand side
	if ~isempty(H1.F)
		H = H1;
		H.F = H1.F\H2.F;
	else
		H = H2;
		H.U12 = H1.A11\H2.U12;
		H.A11 = H1.A11\H2.A11;
		H.U21 = [H2.U21, H1.U21];
		H.V21 = [H2.V21,-H.A11'*H1.V21];
		H.U21 = H1.A22\ H.U21;
		H.A22 = H1.A22\hmatrix_rank_update(H2.A22, -H1.U21 * (H1.V21' * (H1.A11 \ H2.U12) ) ,H2.V12);   
	end
else % case of dense right-hand side
	if ~isempty(H1.F)
		H = H1.F\H2;
	else
		mp = size(H1.A11,2);
		x2 = H1.A11\H2(1:mp,:);
		x1 = H1.A22\(H2(mp+1:end,:) - H1.U21 * (H1.V21' * x2));
		H = [x2; x1];
	end
end
