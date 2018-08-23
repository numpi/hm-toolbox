function B = hss_proper(A)
% Forward stage
if (A.leafnode==1)
    B = A;
    return; % It means that we have called the compression on a HSS matrix with only the root, so nothing to do
end
if (A.A11.leafnode==1)
    [Ut, Pt] = qr(A.A22.U,0);
    [Vt, Qt] = qr(A.A22.V,0);
    A.A22.U = Ut;
    A.A22.V = Vt;
    if (A.topnode==0)
        A.Rr = Pt * A.Rr;
        A.Wr = Qt * A.Wr;
    end
    A.B21 = Pt * A.B21;
    A.B12 = A.B12 * Qt';
    [Ut, Pt] = qr(A.A11.U,0);
    [Vt, Qt] = qr(A.A11.V,0);
    A.A11.U = Ut;
    A.A11.V = Vt;
    if (A.topnode==0)
        A.Rl = Pt * A.Rl;
        A.Wl = Qt * A.Wl;
    end
    A.B21 = A.B21 * Qt';
    A.B12 = Pt * A.B12;
    B = A;
else
    A.A11 = hss_proper(A.A11);
    A.A22 = hss_proper(A.A22);
    
    [Ut, Pt] = qr([A.A11.Rl;A.A11.Rr],0);
    [Vt, Qt] = qr([A.A11.Wl;A.A11.Wr],0);
    j = size(A.A11.Rl,1);
    A.A11.Rl = Ut(1:j,:);
    A.A11.Rr = Ut(j+1:end,:);
    j = size(A.A11.Wl,1);
    A.A11.Wl = Vt(1:j,:);
    A.A11.Wr = Vt(j+1:end,:);
    if (A.topnode==0)
        A.Rl = Pt * A.Rl;
        A.Wl = Qt * A.Wl;
    end
    A.B21 = A.B21 * Qt';
    A.B12 = Pt * A.B12;
    
    [Ut, Pt] = qr([A.A22.Rl;A.A22.Rr],0);
    [Vt, Qt] = qr([A.A22.Wl;A.A22.Wr],0);
    j = size(A.A22.Rl,1);
    A.A22.Rl = Ut(1:j,:);
    A.A22.Rr = Ut(j+1:end,:);
    j = size(A.A22.Wl,1);
    A.A22.Wl = Vt(1:j,:);
    A.A22.Wr = Vt(j+1:end,:);
    if (A.topnode==0)
        A.Rr = Pt * A.Rr;
        A.Wr = Qt * A.Wr;
    end
    A.B21 = Pt * A.B21;
    A.B12 = A.B12 * Qt';
    B = A;
end
end

