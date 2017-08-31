function [ Xu, Xv ] = SylvKrylov(A, B, u, v, k)
%SYLVKRYLOV Summary of this function goes here
%   Detailed explanation goes here
if isa(A, 'hm')
    n = size(A);
    n = n(1);
else
    n = size(A, 1);
end

% k = 10; 

% va = u / norm(u);
% vb = v / norm(v);
[va, RA] = qr(u, 0);
[vb, RB] = qr(v, 0);

% Cs = RA * RB';

QA = va;
QB = vb;

bs = size(va, 2);

while size(QA, 2) < k * bs

    jj = size(QA, 2) / bs;
    
    if jj == 1
        wa = A  * QA(:,1:bs);
        wb = B' * QB(:,1:bs);
    else
        if mod(jj, 2) == 1
            wa = A  * QA(:,end-2*bs+1:end-bs);
            wb = B' * QB(:,end-2*bs+1:end-bs);
        else
            wa = A  \ QA(:,end-2*bs+1:end-bs);
            wb = B' \ QB(:,end-2*bs+1:end-bs);
        end
    end
    
    % Orthogonalization
    wa = wa - QA * (QA' * wa);
    wb = wb - QB * (QB' * wb);
        
    wa = wa - QA * (QA' * wa);
    wb = wb - QB * (QB' * wb);
    
    % norm(wa)
    % norm(wb)
    
    %wa = wa / norm(wa);
    %wb = wb / norm(wb);
    [wa, ~] = qr(wa, 0);
    [wb, ~] = qr(wb, 0);
    
    QA = [ QA , wa ];
    QB = [ QB , wb ];
end

% Theoretical best basis
% [U,~,V] = svd(lyap(A,B,-u*v'));
% QA = U(:,1:k);
% QB = V(:,1:k);

As = QA' * (A * QA);
Bs = QB' * (B * QB);

% Cs = zeros(k); Cs(1,1) = norm(u) * norm(v);
Cs = zeros(k * bs); Cs(1:bs,1:bs) = RA * RB';
% Cs = (QA' * u) * (v' * QB);

Xs = lyap(As, Bs, -Cs);
% reduce solution rank if needed
 [Us,Ss,Vs] = svd(Xs);
 sY = diag(Ss);
 is = sum(abs(sY) / sY(1) > 1e-12);
 Y0 = Us(:,1:is)*diag(sqrt(sY(1:is)));

Xu = -QA * Us(:,1:is) * diag(sqrt(sY(1:is)));
Xv = QB * Vs(:,1:is) * diag(sqrt(sY(1:is)));
% res = norm(A*X + X*B - u*v');
% res = 1;
% X = Xs;

end

