function ind = end(H, i, n)
% Overload end operator for calculating index bounds

if n ~= 2, error('not supported'); end
sz = size(H);
if i == 1,
    ind = sz(1);
elseif i == 2,
    ind = sz(2);
end
