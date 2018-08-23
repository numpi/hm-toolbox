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

