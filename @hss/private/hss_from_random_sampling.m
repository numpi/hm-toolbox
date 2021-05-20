function B = hss_from_random_sampling(obj, Afun, Afunt, Aeval, m, n, target_rank)
%HSS_FROM_RANDOM_SAMPLING Build the HSS representation of a matrix
% 			   using mat-vec multiplication with random (block) vectors
%			   and access to (block) diagonal entries.
%
% B = HSS_RANDOM_SAMPLING is based on the algorithm in [1].
%
% [1] Martinsson, Per-Gunnar. "A fast randomized algorithm for computing a
%     hierarchically semiseparable representation of a matrix." SIAM
%     Journal on Matrix Analysis and Applications 32.4 (2011): 1251-1274

% A few things are not available in Octave, so we use workarounds
% isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

if ~exist('target_rank', 'var')
    target_rank = 0;
end

tol = hssoption('threshold');

B = obj;

if B.topnode == 1 && B.leafnode == 1
    [m, n] = size(B);
    B.D = Aeval((1:m).', 1:n);
    return
end

failed = true;

k = 10;
a = 10; % additional columns used for testing the residual

if target_rank > 0
    k = target_rank;
    a = 0;
end

Ocol = randn(n, k + a);
Scol = Afun(Ocol);
Orow = randn(m, k + a);
Srow = Afunt(Orow);

% Increase on the block size when the required accuracy is not reached
bs = 12;

% Norm estimate through a few iterations of the power method. 
v = randn(n, 1);
converged = false;
nrmold = 0; j = 0;
while ~converged
    v = Afunt(v / norm(v));
    v = Afun(v);
    
    nrm = sqrt(norm(v));
    j = j + 1;
    
    if (nrm - nrmold) < nrmold * 5e-2
        converged = true;
    end
    
    nrmold = nrm;
end    

while failed
    % fprintf('HSS_RANDOM_FROM_SAMPLING :: columns: %d\n', size(Scol, 2));
    [B, failed] = hss_from_random_sampling_rec(B, Aeval, Scol, Srow, ...
        Ocol, Orow, 0, 0, tol, nrm, a);

    if failed
        % fprintf('HSS_FROM_RANDOM_SAMPLING :: Enlarging sampling space to %d\n', k + bs);
        Ocol = [ Ocol, randn(n, bs) ];
        Scol = [ Scol, Afun(Ocol(:, end-bs+1:end))  ];
        Orow = [ Orow, randn(m, bs) ];
        Srow = [ Srow, Afunt(Orow(:,end-bs+1:end)) ];
        k = k + bs;
    end
end

% Construct the rest of the structure
% fprintf('HSS_FROM_RANDOM_SAMPLING :: Recovering the rest of the structure\n');
failed = true;
while failed    
    [B, ~, ~, ~, ~, ~, ~, ~, ~, failed] = ...
        hss_from_random_sampling_rec2(B, Aeval, Scol, Srow, ...
            Ocol, Orow, 0, 0, tol, nrm, a);
        
    if failed
        % fprintf('HSS_FROM_RANDOM_SAMPLING :: Enlarging sampling space to %d\n', k + bs);
        Ocol = [ Ocol, randn(n, bs) ];
        Scol = [ Scol, Afun(Ocol(:, end-bs+1:end))  ];
        Orow = [ Orow, randn(m, bs) ];
        Srow = [ Srow, Afunt(Orow(:,end-bs+1:end)) ];
        k = k + bs;
    end
end

% Remove any temporary data that we might have stored in the leafnodes B12
% anb B21 -- these are not part of the final HSS data structure.
B = clean_structure(B);

end

function [B, failed] = hss_from_random_sampling_rec(B, Aeval, Scol, ...
    Srow, Ocol, Orow, row, col, tol, nrm, a)

if B.leafnode == 1
    failed = false;
    
    [m, n] = size(B.D);
    
    % Only store the matrix D if it has not been initialized yet; in
    % general it might have been initialized in a previous pass.
    if isequal(B.D, zeros(m, n))
        B.D = Aeval(row + 1: row + m, col + 1:col + n);
        B.B12 = [];
        B.B21 = [];
    end
    
    % Compute span of the columns, but only if the old one was not ok
    if isempty(B.B12) || true
        Scol = Scol - B.D * Ocol;
        
        [Q, ~] = colspan(Scol(:, 1:end - a), ...
            sqrt(size(Scol, 2)) * eps, nrm);
        
        Scol2 = Scol(:, end-a+1:end);
        Scol2 = Scol2 - Q * (Q' * Scol2);
        
        %km = 1 + 11 * sqrt(size(Q,2) * size(Scol, 1)); % Martinsson's constant, too pessimistic in practice
        km = 1;
        
        if norm(Scol2) * km > nrm * tol
            failed = true;
            B.B12 = [];
            return;
        end
        
        [Xcol, Jcol] = interpolative(Q', tol);
        
        B.U = Xcol';
        B.B12 = Jcol;
    end
    
    % Compute span of the rows
    if isempty(B.B21) || true
        Srow = Srow - B.D' * Orow;
        [Q, ~] = colspan(Srow(:,1:end-a), ...
            sqrt(size(Srow, 2)) * eps, nrm);
        
        Srow2 = Srow(:,end-a+1:end);
        Srow2 = Srow2 - Q * (Q' * Srow2);
        
        %km = 1 + 11 * sqrt(size(Q,2) * size(Srow, 1)); % Martinsson's constant
        km = 1;
        
        if norm(Srow2) * km > nrm * tol
            failed = true;
            B.B21 = [];
            return;
        end
        
        [Xrow, Jrow] = interpolative(Q', tol);
        
        B.V = Xrow';        
        B.B21 = Jrow;
    end
else
    [B.A11, failed1] = hss_from_random_sampling_rec(B.A11, Aeval, Scol(1:B.ml, :), ...
        Srow(1:B.nl, :), Ocol(1:B.nl, :), Orow(1:B.ml, :), ...
        row, col, tol, nrm, a);
    
    if ~failed1
        [B.A22, failed2] = hss_from_random_sampling_rec(B.A22, Aeval, ...
            Scol(B.ml + 1:end, :), Srow(B.nl + 1:end,:), ...
            Ocol(B.nl + 1:end, :), Orow(B.ml + 1:end, :), ...
            row + B.ml, col + B.nl, tol, nrm, a);
    end
    
    if (failed1 || failed2)
        failed = true;
        return;
    else
        failed = false;
    end
end
end

function [B, Scol, Srow, Ocol, Orow, Jcol, Jrow, U, V, failed] = ...
    hss_from_random_sampling_rec2(B, Aeval, Scol, Srow, Ocol, Orow, row, col, tol, nrm, a)

    failed = false;

    if B.leafnode == 1
        % Recover the information from the previous pass
        U = B.U;
        Jcol = B.B12;
        Scol = Scol(Jcol, :) - B.D(Jcol, :) * Ocol;
        Jcol = row + Jcol;

        V = B.V;
        Jrow = B.B21;
        Srow = Srow(Jrow, :) - B.D(:, Jrow)' * Orow;
        Jrow = col + Jrow;
    else
        [B.A11, Scol1, Srow1, Ocol1, Orow1, Jcol1, Jrow1, U1, V1, failed1] = ...
            hss_from_random_sampling_rec2(B.A11, Aeval, ...
            Scol(1:B.ml, :), Srow(1:B.nl, :), Ocol(1:B.nl, :), ...
            Orow(1:B.ml, :), row, col, tol, nrm, a);
        
        if failed1
            failed = true;
            Scol = []; Srow = []; Ocol = []; Orow = [];
            Jcol = []; Jrow = []; U = []; V = [];
            return;
        end

        [B.A22, Scol2, Srow2, Ocol2, Orow2, Jcol2, Jrow2, U2, V2, failed2] = ...
            hss_from_random_sampling_rec2(B.A22, Aeval, ...
            Scol(B.ml + 1:end, :), Srow(B.nl + 1:end,:), ...
            Ocol(B.nl + 1:end, :), Orow(B.ml + 1:end, :), ...
            row + B.ml, col + B.nl, tol, nrm, a);
        
        if failed2
            failed = true;
            Scol = []; Srow = []; Ocol = []; Orow = [];
            Jcol = []; Jrow = []; U = []; V = [];
            return;
        end

        Ocol2 = V2' * Ocol2;
        Ocol1 = V1' * Ocol1;
        Orow2 = U2' * Orow2;
        Orow1 = U1' * Orow1;

        Jcol = [Jcol1, Jcol2]; Jrow = [Jrow1, Jrow2];
        Ocol = [Ocol1; Ocol2]; Orow = [Orow1; Orow2];

        B.B12 = Aeval(Jcol1, Jrow2);
        B.B21 = Aeval(Jcol2, Jrow1);

        Scol = [Scol1 - B.B12  * Ocol2;  Scol2 - B.B21  * Ocol1 ];
        Srow = [Srow1 - B.B21' * Orow2;  Srow2 - B.B12' * Orow1 ];

        if B.topnode == 0
            % If B.U is non empty we can recover Jcolloc from the previous
            % run, which was successful. The basis is not recomputed.
            if isempty(B.U)
                Q = colspan(Scol(:, 1:end-a), sqrt(size(Scol, 2)) * eps, nrm);

                Scol2 = Scol(:, end-a+1:end);
                Scol2 = Scol2 - Q * (Q' * Scol2);

                %km = 1 + 11 * sqrt(size(Q,2) * size(Scol, 1)); % Martinsson's constant, too pessimistic in practice
                km = 1;

                if norm(Scol2) * km > sqrt(size(Srow, 2)) * nrm * tol
                    failed = true;
                    Scol = []; Srow = []; Ocol = []; Orow = [];
                    Jcol = []; Jrow = []; U = []; V = [];
                    return;
                end

                [Xcol, Jcolloc] = interpolative(Q', tol);
                B.U = Jcolloc;

                B.Rl = Xcol(:, 1:size(Scol1, 1))';
                B.Rr = Xcol(:, size(Scol1, 1)+1:end)';
            else
                Jcolloc = B.U;
            end
            
            Scol = Scol(Jcolloc, :);
            Jcol = Jcol(Jcolloc);
            U = [ B.Rl ; B.Rr ];
            
            % If B.V is non empty we can recover Jrowloc from the previous
            % run, which was successful. The basis is not recomputed.
            if isempty(B.V)
                Q = colspan(Srow(:, 1:end-a), sqrt(size(Srow, 2)) * eps, nrm);

                Srow2 = Srow(:, end-a+1:end);
                Srow2 = Srow2 - Q * (Q' * Srow2);

                %km = 1 + 11 * sqrt(size(Q,2) * size(Scol, 1)); % Martinsson's constant, too pessimistic in practice
                km = 1;

                if norm(Srow2) * km > sqrt(size(Srow, 2)) * nrm * tol
                    failed = true;
                    Scol = []; Srow = []; Ocol = []; Orow = [];
                    Jcol = []; Jrow = []; U = []; V = [];
                    return;
                end

                [Xrow, Jrowloc] = interpolative(Q', tol);
                B.V = Jrowloc;
                B.Wl = Xrow(:, 1:size(Srow1, 1))';
                B.Wr = Xrow(:, size(Srow1, 1)+1:end)';
            else
                Jrowloc = B.V;
            end
            
            Srow = Srow(Jrowloc, :);
            Jrow = Jrow(Jrowloc);
            V = [ B.Wl ; B.Wr ];
        else
            Scol = []; Srow = []; Ocol = []; Orow = [];
            Jcol = []; Jrow = []; U = []; V = [];
        end
    end
end

function [Q, rk] = colspan(S, tol, nrm)
use_qr = false;

if use_qr
    [Q, R, ~] = qr(S, 0);
    rk = sum(abs(diag(R)) > nrm * tol);
    Q = Q(:,1:rk);
else
    [Q, S, ~] = svd(S);
    rk = sum(diag(S) > nrm * tol);
    Q = Q(:,1:rk);
end
end

function y = normest_afun(Afun, Afunt, x, transp)

if strcmp(transp, 'transp')
    y = Afunt(x);
else
    y = Afun(x);
end

end

function B = clean_structure(B)
%CLEAN_STRUCTURE Clean the HSS structure from temporary data.

if B.leafnode
    B.B12 = [];
    B.B21 = [];
else
    B.U = []; B.V = [];
    B.A11 = clean_structure(B.A11);
    B.A22 = clean_structure(B.A22);
end

end
