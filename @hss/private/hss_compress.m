% function that performs the compression of an HSS matrix
function B = hss_compress(A, tol)
% The norm has the side-effect of making the matrix proper
B = hss_proper(A);

% Select the compression strategy according to the user's choice
switch hssoption('compression')
    case 'qr'
        tcomp = @tqr;
    case 'svd'
        tcomp = @tsvd;
    otherwise
        error('unsupported compression method selected');
end

nrmA = norm(A, 2);

if nrmA == 0
    % We cannot request high-relative accuracy if the norm of A is zero
    nrmA = 1;
end

B = backward_stage(B, tol * nrmA, [], [], tcomp);

end

function A = backward_stage(A, tol, S, T, tcomp)
if(A.leafnode==1)
    return
end
if (A.topnode == 1)
    [U,Su,V] = tcomp(A.B12, tol);
    A.B12 = Su; Tl = Su';
    
    if A.A11.leafnode == 0
        A.A11.Rl = A.A11.Rl * U;
        A.A11.Rr = A.A11.Rr * U;
    else
        A.A11.U = A.A11.U * U;
    end
    
    if A.A22.leafnode == 0        
        A.A22.Wl = A.A22.Wl * V;
        A.A22.Wr = A.A22.Wr * V;
    else
        A.A22.V = A.A22.V * V;
    end
    
    [U,Sl,V] = tcomp(A.B21, tol);
    A.B21 = Sl; Tu = Sl';
    
    if A.A11.leafnode == 0
        A.A11.Wl = A.A11.Wl * V;
        A.A11.Wr = A.A11.Wr * V;
        A.A11 = backward_stage(A.A11, tol, Su, Tu, tcomp);
    else
        A.A11.V = A.A11.V * V;
    end
    
    if A.A22.leafnode == 0        
        A.A22.Rl = A.A22.Rl * U;
        A.A22.Rr = A.A22.Rr * U;
        A.A22 = backward_stage(A.A22, tol, Sl, Tl, tcomp);
    else
        A.A22.U = A.A22.U * U;
    end
    
else
    Su = [A.B12, A.Rl * S];
    Tl = [A.B12', A.Wr * T];
    [Us,Su,Vs] = tcomp(Su, tol);
    [Ut,Tl,Vt] = tcomp(Tl, tol);
    k = size(A.B12,2);
    A.B12 = Su * Vs(1:k, :)' * Ut;
    A.Rl = Us' * A.Rl;
    A.Wr = Ut' * A.Wr;
    if A.A11.leafnode == 0
        A.A11.Rl = A.A11.Rl * Us;
        A.A11.Rr = A.A11.Rr * Us;
    else
        A.A11.U = A.A11.U * Us;
	end
    if A.A22.leafnode == 0	
        A.A22.Wl = A.A22.Wl * Ut;
        A.A22.Wr = A.A22.Wr * Ut;
    else
        A.A22.V = A.A22.V * Ut;
    end
    
    Sl = [A.B21, A.Rr * S];
    Tu = [A.B21', A.Wl * T];
    [Us,Sl,Vs] = tcomp(Sl, tol);
    [Ut,Tu,Vt] = tcomp(Tu, tol);
    k = size(A.B21,2);
    A.B21 = Sl * Vs(1:k, :)' * Ut;
    A.Rr = Us' * A.Rr;
    A.Wl = Ut' * A.Wl;
    if A.A22.leafnode == 0
        A.A22.Rl = A.A22.Rl * Us;
        A.A22.Rr = A.A22.Rr * Us;
        A.A22 = backward_stage(A.A22, tol, Sl, Tl, tcomp);
    else        
        A.A22.U = A.A22.U * Us;
	end
	if A.A11.leafnode == 0
        A.A11.Wl = A.A11.Wl * Ut;
        A.A11.Wr = A.A11.Wr * Ut;
        A.A11 = backward_stage(A.A11, tol, Su, Tu, tcomp);
    else        
        A.A11.V = A.A11.V * Ut;
    end
end
end
