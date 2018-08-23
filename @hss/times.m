function C=times(A,B)
if isa(B,'hss')
    if isa(A,'hss')
        C = hss_hadamard_mul(A,B);
    else
        error('Unsupported Hadamard product between hss and dense matrix')
    end
else
    error('Unsupported Hadamard product between hss and dense matrix')
end
end
