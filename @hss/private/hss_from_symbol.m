function H = hss_from_symbol(am, ap, m, n, rowcluster, colcluster)
%HSS_FROM_SYMBOL

if ~exist('m', 'var')
    m = length(am);
    n = length(ap);
end

if am(1) ~= ap(1)
    warning('am(1) ~= ap(1), using ap(1) for the diagonal entries');
end

am = am(:)'; ap = ap(:)';

if (m == n) && (length(am) + length(ap) - 1) < 10 * log(n) ...
        && prod(rowcluster == colcluster) == 1
    bl = length(am) - 1;
    bu = length(ap) - 1;
    
    H = hss('banded', spdiags(ones(m, 1) * [ am(end:-1:2) , ap ], ...
        -bl : bu, m, n), bl, bu, 'cluster', rowcluster, colcluster);
else
    H = hss('handle', ...
        @(v) toepmult_fft(am, ap, m, n, v), ...
        @(v) toepmult_fft(ap, am, n, m, v), ...
        @(I,J) Aeval(am, ap, I, J), m, n, 'cluster', ...
        rowcluster, colcluster);
end

end

function M = Aeval(am, ap, i, j)
if isscalar(i) && isscalar(j)
    M = 0;
    if i > j && i - j < length(am)
        M = am(i - j + 1);
    end
    if j >= i && j - i < length(ap)
        M = ap(j - i + 1);
    end
else
    M = zeros(length(i), length(j));
    
    for ii = 1 : length(i)
        for jj = 1 : length(j)
            M(ii,jj) = Aeval(am, ap, i(ii), j(jj));
        end
    end
end
end

