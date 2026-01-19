function T = telgallery(name, varargin)
%TELGALLERY Generate examples of H-matrices.
%
% Diagonally dominant random matrix ??
% T = TELGALLERY('rand-dd', n)

switch name
    case 'rand-dd'
        T = tel;
        n = varargin{1};
        k = varargin{2};

        % FIXME: We probably want to migrate to 
        % teloption when it is available.
        block_size = hssoption('block-size');        

        if length(varargin) >= 3
            % FIXME: We probably want to be consistent with 
            % the way parameters are passed to @hss, @hodlr.
            c = varargin{3};
        else
            c = build_balanced_cluster(n, block_size);
        end

        % FIXME: Handle degenerate clusters;
        T.D = cell(length(c), 1);
        T.U = cell(length(c), 1);
        T.V = cell(length(c), 1);

        cc = [0, c];
        for j = 1 : length(c)
            sz = cc(j+1) - cc(j);
            T.D{j} = randn(sz, sz);
            T.U{j} = randn(sz, k);
            T.V{j} = randn(sz, k);
        end

        child_cluster = (2 : 2 : length(c)) * k;
        if length(child_cluster) == 1
            T.M = randn(2*k, 2*k);
            T.top = 1;
        else
            T.M = telgallery('rand-dd', ...
                length(c) * k, k, child_cluster);
        end
    otherwise
        error('Unsupported gallery name');
end


end