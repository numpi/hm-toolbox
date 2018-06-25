function rk = qsrank(H)
%QSRANK Obtain the quasiseparability rank of H. 
%
% RK = QSRANK(H) returns the maximum rank in the representation of the
% off-diagonal blocks. 

warning('QSRANK is deprecated, please use HMRANK instead');
rk = hmrank(H);