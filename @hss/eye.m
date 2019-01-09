function I = eye(varargin)
%EYE Create an HSS identity

if isa(varargin{1}, 'hss')
    [rowcluster, colcluster] = cluster(varargin{1});

    I = hss('eye', size(varargin{1}, 1), 'cluster', rowcluster, colcluster);
else
    sz = varargin{1};

    rowcluster = []; colcluster = [];
    if nargin >= 3
        [rowcluster, colcluster] = cluster(varargin{3});
    end

    I = hss('eye', sz(1), 'cluster', rowcluster, colcluster);
end

