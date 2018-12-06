function B = hss_from_random_sampling(Afun, Afunt, Aeval, m, n)
%HSS_FROM_RANDOM_SAMPLING Build the HSS representation of a matrix
% 			   using mat-vec multiplication with random (block) vectors
%			   and access to (block) diagonal entries.
%
% B = HSS_RANDOM_SAMPLING is based on the algorithm in
%
% Martinsson, Per-Gunnar. "A fast randomized algorithm for computing a hierarchically semiseparable representation of a matrix." 
% SIAM Journal on Matrix Analysis and Applications 32.4 (2011): 1251-1274

tol = hssoption('threshold');

B = hss_build_hss_tree(m, n, hssoption('block-size'));

failed = true;

k = 20;

Ocol = randn(n, k);
Scol = Afun(Ocol);
Orow = randn(m, k);
Srow = Afunt(Orow);

while failed
    % fprintf('HSS_RANDOM_FROM_SAMPLING :: columns: %d\n', size(Scol, 2));
    [B, ~, ~, ~, ~, ~, ~, ~, ~, failed] = ...
        hss_from_random_sampling_rec(B, Aeval, Scol, Srow, ...
            Ocol, Orow, 0, 0, tol);
    
    if failed
        Ocol = [ Ocol, randn(n, 10) ];
        Scol = [ Scol, Afun(Ocol(:, end-9:end))  ];
        Orow = [ Orow, randn(m, 10) ];
        Srow = [ Srow, Afunt(Orow(:,end-9:end)) ];
    end
end

B = hss_compress(B, tol);

end

function [B, Scol, Srow, Ocol, Orow, Jcol, ...
          Jrow, U, V, failed] = ...
              hss_from_random_sampling_rec(B, Aeval, Scol, Srow, ...
              Ocol, Orow, row, col, tol)
          
	if B.leafnode == 1
        failed = false;
        
		[m, n] = size(B.D);
        
        % Only store the matrix D if it has not been initialized yet; in
        % general it might have been initialized in a previous pass. 
        if isequal(B.D, zeros(m, n))        
            B.D = Aeval(row + 1: row + m, col + 1:col + n);
        end
        
		Scol = Scol - B.D * Ocol;
		[Q, R, ~] = qr(Scol, 0); rk = sum(abs(diag(R)) > abs(R(1,1)) * eps * size(R, 2));
        
        if rk >= size(Scol, 2) - 10
            failed = true;
        end
        
        Q = Q(:,1:rk);
        
		[Xcol, Jcol] = interpolative(Q', eps);
		B.U = Xcol';
		Scol = Scol(Jcol, :);
        
        U = Xcol';
        Jcol = row + Jcol;

		Srow = Srow - B.D' * Orow;
        
		[Q, R, ~] = qr(Srow, 0); rk = sum(abs(diag(R)) > abs(R(1,1)) * eps * size(R, 2));
                
        if rk >= size(Srow, 2) - 10
            failed = true;
        end
        
        Q = Q(:,1:rk);
        
		[Xrow, Jrow] = interpolative(Q.', eps);

		B.V = Xrow.';
		Srow = Srow(Jrow, :);
        
        V = Xrow.';
        
        Jrow = col + Jrow;
	else
		[B.A11, Scol1, Srow1, Ocol1, Orow1, Jcol1, Jrow1, U1, V1, failed1]  = hss_from_random_sampling_rec(B.A11, Aeval, Scol(1:B.ml, :), Srow(1:B.nl, :), Ocol(1:B.nl, :), Orow(1:B.ml, :), row, col, tol);
		[B.A22, Scol2, Srow2, Ocol2, Orow2, Jcol2, Jrow2, U2, V2, failed2]  = hss_from_random_sampling_rec(B.A22, Aeval, Scol(B.ml + 1:end, :), Srow(B.nl + 1:end,:), Ocol(B.nl + 1:end, :), Orow(B.ml + 1:end, :), row + B.ml, col + B.nl, tol);
        
        if (failed1 || failed2)
            failed = true;
            
            Scol = []; Srow = []; Ocol = []; Orow = []; 
            Jcol = []; Jrow = []; U = []; V = [];
            
            return;
        else
            failed = false;
        end
        
        B.B12 = Aeval(Jcol1, Jrow2);
		B.B21 = Aeval(Jcol2, Jrow1);
        
        if B.topnode == 0
            Ocol2 = V2' * Ocol2;
            Ocol1 = V1' * Ocol1;
            Orow2 = U2' * Orow2;
            Orow1 = U1' * Orow1;

            Jcol = [Jcol1, Jcol2]; Jrow = [Jrow1, Jrow2];
            Scol = [Scol1 - B.B12  * Ocol2;  Scol2 - B.B21  * Ocol1 ]; 
            Srow = [Srow1 - B.B21' * Orow2;  Srow2 - B.B12' * Orow1 ];
            Ocol = [Ocol1; Ocol2]; Orow = [Orow1; Orow2];

            [Q, R, ~] = qr(Scol, 0); rk = sum(abs(diag(R)) > abs(R(1,1)) * eps * size(R, 2));

            Q = Q(:,1:rk);
            
            [Xcol, Jcolloc] = interpolative(Q', eps);

            B.Rl = Xcol(:, 1:size(Scol1, 1))';
            B.Rr = Xcol(:, size(Scol1, 1)+1:end)';
            Scol = Scol(Jcolloc, :);
            Jcol = Jcol(Jcolloc);
            U = [ B.Rl ; B.Rr ];
            
            [Q, R, ~] = qr(Srow, 0); rk = sum(abs(diag(R)) > abs(R(1,1)) * eps * size(R, 2));

            Q = Q(:,1:rk);
            
            [Xrow, Jrowloc] = interpolative(Q', eps);

            B.Wl = Xrow(:, 1:size(Srow1, 1))';
            B.Wr = Xrow(:, size(Srow1, 1)+1:end)';
            Srow = Srow(Jrowloc, :);
            Jrow = Jrow(Jrowloc);
            V = [ B.Wl ; B.Wr ];            
        else
            Scol = []; Srow = []; Ocol = []; Orow = []; 
            Jcol = []; Jrow = []; U = []; V = [];
        end
	end
end

