function c = build_balanced_cluster(n, block_size)

if n <= block_size
    c = n;
    return;
else
    [n1, n2] = split_indices(n);

    left  = build_balanced_cluster(n1, block_size);
    right = n1 + build_balanced_cluster(n2, block_size);

    if length(left) < length(right)
        missing = length(right) - length(left);
        left = [ left, left(end) * ones(1, missing) ];
    end

    if length(right) < length(left)
        missing = length(left) - length(right);
        right = [ right(1) * ones(1, missing), right ];
    end

    c = [left, right];
end

end

function [n1, n2] = split_indices(n)
    n1 = ceil(n / 2);
    n2 = n - n1;
end