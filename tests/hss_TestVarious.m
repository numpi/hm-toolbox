function hss_TestVarious

hssoption('block-size', 256);
hssoption('threshold', 1e-12);
tol = hssoption('threshold');

% Test the Toeplitz solver
for j = 1 : 2
    n = 16384;
    
    c = rand(n, 1);
    r = rand(n, 1).';
    
    if j == 1
        realcomplex = 'real';
    else
        c = c + 1i * rand(n, 1);
        r = r + 1i * rand(n, 1).';
        realcomplex = 'complex';
    end
    
    c(1) = r(1);
    scl = norm(c) + norm(r);
    c = c / scl;
    r = r / scl;
    
    v = rand(n, 1);
    if j == 2
        v = v + 1i * rand(n, 1);
    end
    
    tic;
    x = toeplitz_solve(c, r, v);
    t = toc;
    
    % Try the dense solver -- that should be slower!
    tic; toeplitz(c, r) \ v; td = toc;
    
    CheckTestResult(t, '<', td, ...
        sprintf('Performances of the Toeplitz fast solver (%s case)', realcomplex));
    
    % Check the accuracy -- here we let it a little bit loose on the tolerance,
    % because at the moment the Martinsson's constructor does not guarantee the
    % relative accuracy with tol. We estimate the norm of the Toeplitz part by
    % the Frobenius norm for simplicity.
    nrm_toep = sqrt(sum((n:-1:1) .* abs(c.').^2) + sum((n-1:-1:1) .* abs(r(2:end)).^2));
    CheckTestResult(norm(toeplitz(c, r) * x - v) / (norm(v) + nrm_toep * norm(x)), '<', ...
        1e3 * sqrt(n) * tol, ...
        sprintf('Accuracy of the Toeplitz solver (%s case)', realcomplex));
end


end

function y = toepmult_afun(am, ap, m, n, x, trasp)
%

if strcmp(trasp, 'trasp')
    tmp = am;
    am = ap;
    ap = tmp;
    tmp = n;
    n = m;
    m = tmp;
end

y = toepmult_fft(am, ap, m, n, x);

end

function u = toepmult_fft(am, ap, m, n, v)
%TOEPMULT_FFT Fast multiplication of a Toeplitz matrix times a vector.
%
% U = TOEPMULT_FFT(AM, AP, M, N, V) computes the matrix-matrix product A*V
%     where A is a Toeplitz matrix defined by the vectors AM and AP and has
%     size M x N, where N = length(V).
%

realflag = isreal(am) && isreal(ap) && isreal(v);

n1=length(am);
n2=length(ap);

if (size(v, 1) == n)
    mn1 = min(n1,m);
    mn2 = min(n,n2);
    
    N = max(m + mn2, n + mn1);
    a = zeros(1,N);
    a(1:mn1) = am(1:mn1);
    a(end:-1:end-mn2+2) = ap(2:mn2);
    
    w = zeros(N,size(v,2));
    w(1:n,:) = v;
    wf=fft(w);
    af=fft(a);
    
    u = zeros(length(af), size(wf, 2));
    af = reshape(af, length(af), 1);
    
    % Choose the fastest possible way to perform the multiplication
    if size(wf, 1) > size(wf, 2)
        for i = 1 : size(wf, 2)
            u(:,i) = af .* wf(:,i);
        end
    else
        for i = 1 : size(wf, 1)
            u(i,:) = af(i) * wf(i,:);
        end
    end
    
    u = ifft(u);
    u = u(1:m,:);
    
    if (realflag)
        u = real(u);
    end
else
    error('Size mismatch in matrix vector multiplication');
end

end

