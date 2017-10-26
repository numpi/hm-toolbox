function [ Xu, Xv ] = SylvKrylov(A, B, u, v, k)
%SYLVKRYLOV Summary of this function goes here
%   Detailed explanation goes here

n = size(A);
n = n(1);

bs = size(v, 2);

% If the problem is small enough, solve it using the dense solver
if n <= bs * k
	%fprintf('dense\n')
	[Xu, Xs, Xv] = svd(lyap(full(A), full(B), u*v'));
	sY = diag(Xs); is = sum(abs(sY) / sY(1) > 1e-12); Xu = Xu(:,1:is) * diag(sqrt(sY(1:is))); Xv = Xv(:,1:is)* diag(sqrt(sY(1:is)));
	return;
end
%fprintf('NON dense\n')
% k = 10; 

% va = u / norm(u);
% vb = v / norm(v);
[va, RA] = qr(u, 0);
[vb, RB] = qr(v, 0);


% Cs = RA * RB';

QA = va;
QB = vb;


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
    
    % norm(wb)
    
    %wa = wa / norm(wa);
    %wb = wb / norm(wb);
    [wa, ~] = qr(wa, 0);
    [wb, ~] = qr(wb, 0);
    
    %[QA,~] = qr([ QA , wa ], 0);
    %[QB,~] = qr([ QB , wb ], 0);
	QA = [ QA, wa ];
	QB = [ QB, wb ];
end

% Theoretical best basis
% [U,~,V] = svd(lyap(A,B,-u*v'));
% QA = U(:,1:k);
% QB = V(:,1:k);

[QA,~] = qr(QA, 0);
[QB,~] = qr(QB, 0);

% norm(QA' * QA - eye(size(QA,2)))

As = QA' * (A * QA);
Bs = QB' * (B * QB);

% Cs = zeros(k); Cs(1,1) = norm(u) * norm(v);
% Cs = zeros(k * bs); Cs(1:bs,1:bs) = RA * RB';
Cs = (QA' * u) * (QB' * v)';
% Cs = (QA' * u) * (v' * QB);

Xs = lyap(As, Bs, -Cs);
% [Us,Ss,Vs] = svd(Xs);
% reduce solution rank if needed 
[Us,Ss,Vs] = svd(Xs); sY = diag(Ss); is = sum(abs(sY) / sY(1) > 1e-12); Y0 = Us(:,1:is)*diag(sqrt(sY(1:is)));
Xu = -QA * Us(:,1:is) * diag(sqrt(sY(1:is))); Xv = QB * Vs(:,1:is) * diag(sqrt(sY(1:is)));

% res = norm(A*X + X*B - u*v');
% res = 1;
% X = Xs;

end