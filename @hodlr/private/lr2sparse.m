function A = lr2sparse(U, V, tol)
% Sparse (approximation) of factorized matrix U*V'
    
m = size(U,1); n = size(V,1);
if m*n == 0, return, end

if tol == 0,
    % Only consider submatrix of structurally zero entries
    A = sparse(m,n);
    i = find( any(U,2) ); j = find( any(V,2) );
    A(i,j) = sparse( U(i,:)*V(j,:)' );
else
% Optional step 0: Equilibrate U and V (sometimes helps)
%    [U,RU] = qr(U,0);
%    [V,RV] = qr(V,0);
%    [QU, D, QV] = svd( RU*RV');
%    D = sqrt(D);
%    U = U*QU*D;
%    V = V*QV*D;
    
    % Step 1: Find entries (i,j) for which
    % norm(U(i,:))*norm(V(i,:)) > tol.
    
    normsU = sqrt(sum(U.^2,2));
    normsV = sqrt(sum(V.^2,2));

    [normsU,permi] = sort(normsU, 'descend');
    [normsV,permj] = sort(normsV, 'descend');

    i = []; j = [];
    for col = 1:size(V,1),
        % Find first entry <= tol / normsV(col)
        ind = find( normsU <= tol / normsV(col), 1 );
        if isempty( ind ), ind = length( normsU ) + 1; end
        ind = ind - 1;
        j = [ j; col*ones(ind,1) ];
        i = [ i; (1:ind)' ];
        % Shorten vector
        normsU = normsU(1:ind);
        if ind==0, break, end
    end
    
    % Resort indices
    i = permi(i);
    j = permj(j);
    [~,ind] = sort(j+1i*i,'ComparisonMethod','real');
    i = i(ind);
    j = j(ind);
    
    % Step 2: Evaluate all remaining entries column by column
    aj = accumarray(j,1);
    vals = zeros(length(i),1);
   
    ind = 1;
    for entries = aj',
        % current row indices and column index
        row = i(ind:ind+entries-1); col = j(ind);
        vals(ind:ind+entries-1) = U(row,:) * V(col,:)';
        ind = ind + entries;
    end
    
    % Step 3: Treshold remaining entries and create sparse matrix
    ind = find( abs(vals)>tol );
    vals = vals(ind);
    i = i(ind); j = j(ind);
    
    A = sparse( i, j, vals, m, n );    
end
