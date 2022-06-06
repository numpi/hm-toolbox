function A = hss_force_symmetry(A, ud)

	A = proper(A);
	
	A = hss_force_symmetry_ric(A, ud);
end

function A = hss_force_symmetry_ric(A, ud)	
	
	if A.leafnode == 1
		if strcmp(ud, 'down')
			A.D = tril(A.D) + tril(A.D, -1)';
		else
			A.D = triu(A.D) + triu(A.D, 1)';
		end
		return
	end

	if strcmp(ud, 'down')
		if A.A11.leafnode == 1
			A.A11.D = tril(A.A11.D) + tril(A.A11.D, -1)';
			A.A22.D = tril(A.A22.D) + tril(A.A22.D, -1)';

			X21 = A.A11.V' * A.A11.U;
		%   X12 = A.A22.V' * A.A22.U;
			
			A.A11.V = A.A11.U;
			A.A22.V = A.A22.U;
			
		%	A.B12 = A.B12 * X12;
			A.B21 = A.B21 * X21;
			
			if A.topnode ~= 1
				A.Wl = A.Rl;
				A.Wr = A.Rr;
			end
			
			A.B12 = A.B21';

		else
			A.A11 = hss_force_symmetry_ric(A.A11, ud);
			A.A22 = hss_force_symmetry_ric(A.A22, ud);
			
			if A.topnode ~= 1
				A.Wl = A.Rl;
				A.Wr = A.Rr;
			end
			
			A.B12 = A.B21';
		end
	else
		if A.A11.leafnode == 1
			A.A11.D = triu(A.A11.D) + triu(A.A11.D, 1)';
			A.A22.D = triu(A.A22.D) + triu(A.A22.D, 1)';

		    X12 = A.A11.V' * A.A11.U;
			
			A.A11.V = A.A11.U;
			A.A22.V = A.A22.U;
			
			A.B12 = X12 * A.B12;
			
			if A.topnode ~= 1
				A.Rl = A.Wl;
				A.Rr = A.Wr;
			end
			
			A.B21 = A.B12';

		else
			A.A11 = hss_force_symmetry_ric(A.A11, ud);
			A.A22 = hss_force_symmetry_ric(A.A22, ud);
			
			if A.topnode ~= 1
				A.Rl = A.Wl;
				A.Rr = A.Wr;
			end
			
			A.B21 = A.B12';
		end	
	
	
	end	
end
