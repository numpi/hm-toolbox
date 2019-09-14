function nrm = normest_Afun(Afun, Afunt, n)
s = 0;
v = randn(n, 1);
for i = 1 : 10
    olds = s;
    s = norm(v);
    if abs(olds - s) < abs(s) * 1e-3 || s == 0
        break;
    end
    
    v = v / s;
    w = v;
    w = Afun(w);
    w = Afunt(w);
    v = w;
end

nrm = sqrt(s);
end
