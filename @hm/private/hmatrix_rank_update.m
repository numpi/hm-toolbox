function H = hmatrix_rank_update(H, U, V)
%HMATRIX_RANK_UPDATE Perform a low rank update to H. 

H = H + hm('low-rank', U, V);

