function T = telgallery(name, varargin)
%TELGALLERY Generate examples of H-matrices.
%
%
% T = TELGALLERY('rand-dd',n,k,rho) generates a n x n matrix with the
% following telescopic decomposition:
%D^(L-i) is diagonal with spectrum in [rho^(2i+1),rho^(2i)] with s<1
%U^(L-i), V^(L-i) rectangular blockdiag matrices with random normal
%(column-wise) blocks
switch name
    case 'rand-dd'
        if isempty(varargin)
            n = 2048; k = 10; rho = 0.5;
        elseif length(varargin) == 1
            n = varargin{1}; k = 10; rho = 0.5;
        elseif length(varargin) == 2
            n = varargin{1}; k = varargin{2}; rho = 0.5;
        elseif length(varargin) == 3
            n = varargin{1}; k = varargin{2}; rho = varargin{3};
        end
        T = tel;
        
        % FIXME: We probably want to migrate to 
        % teloption when it is available.
        block_size = hssoption('block-size');        

        if length(varargin) >= 4
            % FIXME: We probably want to be consistent with 
            % the way parameters are passed to @hss, @hodlr.
            c = varargin{4};
        else
            c = build_balanced_cluster(n, block_size);
        end
        % FIXME: Handle degenerate clusters;
        T.D = cell(length(c), 1);
        T.U = cell(length(c), 1);
        T.V = cell(length(c), 1);

        cc = [0, c];
        l = log2(length(c)); %FIXME: handle degenerate clusters depth
        for j = 1 : length(c)
            sz = cc(j+1) - cc(j);
            d = rho^(2*l + 1) + (rho^(2*l) - rho^(2*l + 1)) * rand(sz,1);
            D = diag(d);
            T.D{j} = D;
            U = randn(sz,k);  
            [U,~] = qr(U,0);  
            T.U{j} = U;
            V = randn(sz,k);  
            [V,~] = qr(V,0); 
            T.V{j} = V;
        end
        child_cluster = (2 : 2 : length(c)) * k;
        if length(child_cluster) == 1
            d = rho^(4*l + 1) + (rho^(4*l) - rho^(4*l + 1)) * rand(2*k,1);
            T.M = diag(d);
            T.top = 1;
        else
            T.top = 0;
            T.M = telgallery('rand-dd', ...
                length(c) * k, k, child_cluster);
        end
    otherwise
        error('Unsupported gallery name');
end


end