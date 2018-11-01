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
