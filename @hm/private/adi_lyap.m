function X = adi_lyap(A, C, tol, maxit)
%ADI_LYAP Solve the Lyapunov equation using the ADI

n = size(A,1);
I = hm('diagonal', ones(n,1));

% compute the ADI shifts to obtain a residual norm <tol
if norm(A-A','fro') < 1e-12
    alpha = 0;
    %opts.tol = 1e-5;
    a = eigs(-A, 1, 'opts','sm' );
    b = eigs(-A, 1, 'opts','lm' );
else
    ll = eig(-full(A));
    a = min(abs(ll));
    b = max(abs(ll));
    alpha = max(atan(abs(imag(ll))./abs(real(ll))));
end

% Compute the number of shifts that guarantees to achieve the tolerance tol
p = myadipars(a,b,alpha,tol);
% otherwise fix the number of shifts (good results can be obtained anyway)
%nshift=input('How many shifts do you want to compute?\n');
%p=myadipars_n(a,b,alpha,nshift);

J = length(p);
%disp(['Number of computed shifts: ', num2str(J)])
j = 1;

% zero initial guess
X = hm('diagonal', zeros(n,1));

for i = 1:maxit
    
    % compute the current sparse approximate inverse
    AshiftedInv = (A+p(j)*I)\I;
    Ashifted = A-p(j)*I;
    
    X = AshiftedInv * C - AshiftedInv * X * Ashifted;
    X = AshiftedInv * C - AshiftedInv * X' * Ashifted;
    
    res=norm(A * X + X * A - C,'fro') / norm(C,'fro');
    %fprintf('It: %d, norm_R: %10.5e\n',i, res)
    %fprintf('*************************\n')
    
    if res<tol
        break
    end
    
    % re use the shifts
    if j == J
        j = 1;
    else
        j = j + 1;
    end
    
    %pause
end
X = -X;
