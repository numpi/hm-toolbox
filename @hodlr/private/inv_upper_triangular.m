function IH = inv_upper_triangular(H)
IH = H;
IH.A11 = inv(H.A11);
IH.A22 = inv(H.A22);
IH.U12 = hodlr_mtimes_dense(-IH.A11 , H.U12);
IH.V12 = hodlr_mtimes_dense(IH.A22' , H.V12);
