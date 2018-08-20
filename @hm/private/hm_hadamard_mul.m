function C = hm_hadamard_mul(A, B);
tol = hssoption('threshold');
C = hm_hadamard_mul_ric(A, B);
C = compress_hmatrix(C, tol);
end

function C = hm_hadamard_mul_ric(A, B)
C = hm();
	if xor(isempty(A.F), isempty(B.F))
		error('HM_HADAMARD_MUL:: the two hm matrices have not compatible partitioning')
	end
	if ~isempty(A.F)
		C.F = A.F .* B.F;
	else
		if size(A.U21, 1) ~= size(B.U21, 1) || size(A.U12, 1) ~= size(B.U12, 1)|| ...
			size(A.V12, 1) ~= size(B.V12, 1) || size(A.V21, 1) ~= size(B.V21, 1)
			error('HM_HADAMARD_MUL:: the two hm matrices have not compatible partitioning')
		end
		for j = 1:size(A.U12,2) 
			C.U12 = [C.U12, diag(A.U12(:,j)) * B.U12];
            		C.V12 = [C.V12, diag(A.V12(:,j)) * B.V12];			
		end
		for j = 1:size(A.U21,2)  
			C.U21 = [C.U21, diag(A.U21(:,j)) * B.U21];
			C.V21 = [C.V21, diag(A.V21(:,j)) * B.V21];
		end
		C.A11 = hm_hadamard_mul_ric(A.A11, B.A11);
		C.A22 = hm_hadamard_mul_ric(A.A22, B.A22);
	end
end
