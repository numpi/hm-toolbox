function I = eye(varargin)
%EYE Create an HODLR identity

if isa(varargin{1}, 'hodlr')
    [rowcluster, colcluster] = cluster(varargin{1});

    I = hodlr('eye', size(varargin{1}, 1), 'cluster', rowcluster, colcluster);
else
    sz = varargin{1};

    rowcluster = []; colcluster = [];
    % We copy the cluster only if the third argument has the right
	% dimensions -- otherwise we use it just as an indication of the type.
    if nargin >= 3 && all(sz == size(varargin{3}))
        [rowcluster, colcluster] = cluster(varargin{3});
    end

    I = hodlr('eye', sz(1), 'cluster', rowcluster, colcluster);
end

