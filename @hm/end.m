function ind = end(H,i,n)
% Overload end operator for calculating index bounds

if n~=2, error('not supported'); end
if i==1,
   ind = H.sz(1);
elseif i==2,
   ind = H.sz(2);
end
