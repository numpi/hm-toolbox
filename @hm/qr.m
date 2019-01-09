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

if m~=n | ~check_cluster_equality(A), error('Input matrix must be square and have a square partition.'); end
BL = zeros(0,0); BR = zeros(0,n);
C = zeros(0,n);    
Anrm = norm(A);

[Y, YBL, YBR, YC, T, A] = qr_iter(A, BL, BR, C, Anrm);

if nargout <= 2,
    [r, c] = cluster(Y);
    Q = hmatrix_minus( hm('eye', n, 'cluster', r), Y*T*Y');
    Y = Q;
    T = A;
end

end

function [YA, BL, YBR, YC, T, A] = qr_iter(A, BL, BR, C, Anrm)

m = size(A,1);
q = size(C,1);

if size(BL,1)>0,
   [BL,R] = qr(BL,0);
   BR = R*BR;
end
p = size(BR,1);


% !!! The following check should be replaced by a proper check for leaf.
if ~isempty(A.F),
    [Y, T, R] = qrWY([A.F;BR;C]);
    
    YA  = hm( Y(1:m,:), 'cluster', m );
    YBR = Y(m+1:m+p,:);
    YC  = Y(m+p+1:end,:);
    A = hm( R(1:m,:), 'cluster', m );
    T = hm( T );
else
    % Compute QR decomposition of first block column
    [m1,n1] = size(A.A11);
    BC = [BR;C];
    [YA11, YBL1, YBR1, YC1, T1, A.A11] = qr_iter(A.A11, A.U21, A.V21', BC(:,1:n1), Anrm);
    SL = [YA11'*A.U12, YBR1', YC1'];
    SR = [A.V12, A.A22'*YBL1, BC(:,n1+1:end)'];
    [ SL, SR ] = compress_factors(SL, SR, Anrm);
    SL = T1'*SL;

    % Update second block column
    [ A.U12, A.V12 ] = compress_factors( [A.U12, -YA11*SL] , [A.V12, SR], Anrm );
    A.A22 = hmatrix_rank_update(A.A22, YBL1, -SR* ( YBR1 * SL )', Anrm);
    BC(:,n1+1:end) = BC(:,n1+1:end) - (YC1*SL)*SR';

    % Compute QR decomposition of second block column
    [m2,n2] = size(A.A22);
    
    [YA22, YBL2, YBR2, YC2, T2, A.A22] = qr_iter(A.A22, zeros(0,0), zeros(0,n2), BC(:,n1+1:end), Anrm);
    
    % Set Y
    YA = hm; YA.A11 = YA11; YA.A22 = YA22; YA.sz = size(A);
    YA.U21 = YBL1; YA.V21 = YBR1';
    YA.U12 = zeros(m1,0); YA.V12 = zeros(n2,0);
    
    YBR = [ YC1(1:p,:), YC2(1:p,:) ];
    YC = [ YC1(p+1:end,:), YC2(p+1:end,:) ];
    
    % Clean up A
    A.U21 = zeros(m2,0); A.V21 = zeros(n1,0);
    
    % Set T
    T = hm; T.sz = [n1+n2,n1+n2];
    T.A11 = T1; T.A22 = T2;
    T.U21 = zeros(n2,0); T.V21 = zeros(n1,0);
    T12L = [ YBR1', YC1(1:p,:)', YC1(p+1:end,:)' ];
    T12R = [ YA22'*YBL1, YBR2', YC2' ];
    [ T12L, T12R ] = compress_factors( T12L, T12R, 1 );
    T.U12 = -T1*T12L;
    T.V12 = T2'*T12R; 
end   

end
