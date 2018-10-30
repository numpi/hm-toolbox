function F = chol(A)
%ULV   Computes the ULV factorization of A.
%
%      F = CHOL(A) returns a structure containing a
%      parametrization of the CHOL of A.
%      This can be used to compute X = A\B with the command
%      X = CHOL_SOLVE(F, B);
if ~issymmetric(A)
	A = hss_symmetrize(A);
end
F = hss_chol_fact(A);
end
