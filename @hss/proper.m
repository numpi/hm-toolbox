function H = proper(H)
% Return a proper hss representation of the matrix H, i.e. one where all the U,V,W and R factors are orthogonal
	H = hss_proper(H);
end 
