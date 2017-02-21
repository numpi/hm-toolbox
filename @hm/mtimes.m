function H = mtimes(H1, H2)
%MTIMES Matrix multiplication

% Multiplication H * v
if isfloat(H2)
	if isscalar(H2)
		if ~isempty(H1.F)
			H = hm();
			H.F = H1.F * H2;
			H.sz = H1.sz;
		else
			H = H1;
			H.A11 = H1.A11 * H2;
			H.A22 = H1.A22 * H2;
			H.U21 = H1.U21 * H2;
			H.U12 = H1.U12 * H2;
		end
	else
		if ~isempty(H1.F)
			H = H1.F * H2;
		else
			mp = H1.A11.sz(2);
			Hu = H1.A11 * H2(1:mp,:) + H1.U12 * (H1.V12' * H2(mp+1:end,:));
			Hl = H1.U21 * (H1.V21' * H2(1:mp,:)) + H1.A22 * H2(mp+1:end,:);		
			H = [ Hu ; Hl ];
		end
	end
	
	return;
end

% Multiplication w' * H
if isfloat(H1)
	if ~isempty(H2.F)
		H = H1 * H2.F;
	else
		if isscalar(H1)
			H = H2 * H1;
		else
			mp = H2.A11.sz(1);
			Hl = H1(:,1:mp) * H2.A11 + (H1(:,mp+1:end) * H2.U21) * H2.V21';
			Hr = (H1(:,1:mp) * H2.U12) * H2.V12' + H1(:,mp+1:end) * H2.A22;
			H = [ Hl , Hr ];
		end
	end
	
	return;
end

% Multiplication of two H-matrices
H = H1;
if ~isempty(H1.F)
	H.F = H1.F * H2.F;
else
	H = H1;
	
	% H.A11 = H1.A11 * H2.A11 + hm('low-rank', H1.U12 * (H1.V12' * H2.U21), H2.V21);
	% H.A22 = hm('low-rank', H1.U21 * (H1.V21' * H2.U12), H2.V12) + H1.A22 * H2.A22;
		
	H.A11 = hmatrix_rank_update(H1.A11 * H2.A11, H1.U12 * (H1.V12' * H2.U21), H2.V21);
	H.A22 = hmatrix_rank_update(H1.A22 * H2.A22,  H1.U21 * (H1.V21' * H2.U12), H2.V12);
	
	[H.U12, H.V12] = compress_factors([ H1.A11 * H2.U12, H1.U12 ], ...
									 [ H2.V12, (H1.V12' * H2.A22)']);
	[H.U21, H.V21] = compress_factors([ H1.A22 * H2.U21, H1.U21 ], ...
									 [ H2.V21, (H1.V21' * H2.A11)']);								 
end

end

