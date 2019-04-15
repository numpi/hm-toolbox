function I = eye(varargin)
%EYE Create an HM identity

if isa(varargin{1}, 'hm')
    [rowcluster, colcluster] = cluster(varargin{1});

    I = hm('eye', size(varargin{1}, 1), 'cluster', rowcluster, colcluster);
else
    sz = varargin{1};

    rowcluster = []; colcluster = [];
    % We copy the cluster only if the third argument has the right
	% dimensions -- otherwise we use it just as an indication of the type.
    if nargin >= 3 && all(sz == size(varargin{3}))
        [rowcluster, colcluster] = cluster(varargin{3});
    end

    I = hm('eye', sz(1), 'cluster', rowcluster, colcluster);
end

