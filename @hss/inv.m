function B = inv(A)
%INV Compute the inverse of an HSS matrix.

if size(A, 1) ~= size(A, 2)
    error('This matrix is not square');
end

B = A \ eye(size(A), 'like', A);

end
