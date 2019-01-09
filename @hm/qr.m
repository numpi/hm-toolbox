function [Y, T, A] = qr(A)
%QR Computes QR decomposition of HODLR matrix.
%
% [Q,R] = QR(A) computes a QR decomposition A = Q*R of a HODLR matrix
%     such that Q is a numerically orthogonal HODLR matrix and R is an upper
%     triangular matrix.
%
% [Y,T,R] = QR(A) returns the factor Q in terms of its compact WY
%     representation: Q = I-Y*T*Y' with an upper triangular matrix T and a
%     lower triangular matrix Y.
%
% Currently, only square matrices are supported.
%
% [1] D. Kressner and A. Šusnjara. (2018). Fast QR decomposition of HODLR
%     matrices. Technical report, September 2018.

[m,n] = size(A);
if m~=n, error('Input matrix must be square'); end
BL = zeros(0,0); BR = zeros(0,n);
C = zeros(0,n);    
[Y, YBL, YBR, YC, T, A] = qr_iter(A, BL, BR, C)

end

function [YA, BL, YBR, YC, T, A] = qr_iter(A, BL, BR, C)

m = size(A,1);
p = size(BR,1);
q = size(C,1);

if size(BL,1)>0,
   [BL,R] = qr(BL,0);
   BR = R*BR;
end

% !!! The following check should be replaced by a proper check for leaf.
if ~isempty(A.F),
    [Y, T, R] = qrWY([A.F;BR;C]);
    
    YA  = Y(1:m,:);
    YBR = Y(m+1:m+p,:);
    YC  = Y(m+p+1:end,:);
    A = R(1:m,:);
else
    % Compute QR decomposition of first block column
    n1 = size(A.A11, 2);
    size(BR(:,1:n1))
    size(C(:,1:n1))
    [Y1, T1, R1] = qr_iter(A.A11, A.U21, A.V21', [ BR(:,1:n1); C(:,1:n1) ]);
    Y1
    T1
    
disp('hallo')
end
    

end
