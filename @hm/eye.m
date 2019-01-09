function I = eye(varargin)
%EYE Create an HM identity

if isa(varargin{1}, 'hm')
    I = hm('eye', size(varargin{1}, 1), 'cluster', cluster(varargin{1}));
else
    sz = varargin{1};

    hcluster = [];
    if nargin >= 3
        hcluster = cluster(varargin{3});
    end

    I = hm('eye', sz(1), 'cluster', hcluster);
end

