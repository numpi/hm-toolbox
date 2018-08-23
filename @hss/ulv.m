function F = ulv(A)
%ULV   Computes the ULV factorization of A.
%
%      F = ULV(A) returns a structure containing a
%      parametrization of the ULV of A.
%      This can be used to compute X = A\B with the command
%      X = ULV_SOLVE(F, B);
F = hss_ulv_fact(A);
end
