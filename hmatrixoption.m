function opt = hmatrixoption(key, value)
%HMATRIXOPTION Set or get an option for the hmatrix toolbox.
%
% Valid options are:
%   'block-size': Integer representing the minimum block size.
%   'threshold': Value used for off-diagonal truncation.
%   'compression': String indicating the strategy of compression (qr or svd)
%   'norm': Norm to use in truncations (2 or 'fro')

global hmatrix_block_size
global hmatrix_threshold
global hmatrix_compression
global hmatrix_norm

% Private options -- not exposed to users

% This option sets the maximum allowed rank for low-rank blocks with
% respect the dimension: for a block of size (m,n) the maximum allowed rank
% is hmatrixoption('rank-ratio') * min(m,n)
global hmatrix_rank_ratio

if isempty(hmatrix_rank_ratio)
    hmatrix_rank_ratio = 0.5;
end

if isempty(hmatrix_block_size)
    hmatrix_block_size = 256;
end

if isempty(hmatrix_compression)
    hmatrix_compression = 'svd';
end

if isempty(hmatrix_threshold)
    hmatrix_threshold = 1e-12;
end

if isempty(hmatrix_norm)
    hmatrix_norm = 2;
end

if ~exist('key', 'var')
    error('Please specify a key');
end

if ~exist('value', 'var')
    switch key
        case 'block-size'
            opt = hmatrix_block_size;
        case 'threshold'
            opt = hmatrix_threshold;
        case 'compression'
            opt = hmatrix_compression;
        case 'norm'
            opt = hmatrix_norm;
        case 'rank-ratio'
            opt = hmatrix_rank_ratio;
        otherwise
            error('Unsupported option specified');
    end
else
    switch key
        case 'block-size'
            if value <= 2
                error('minimum block size must be at least 3');
            else
                hmatrix_block_size = value;
            end
        case 'threshold'
            if value < 0
                error('threshold has to be positive');
            else
                hmatrix_threshold = max(eps, value);
            end
        case 'norm'
            if (isscalar(value) && value ~= 2) && ~strcmp(value, 'fro')
                error('Unsupported norm specified');
            else
                hmatrix_norm = value;
            end
        case 'compression'
            if strcmp(value,'qr') || strcmp(value, 'svd')
                hmatrix_compression = value;
            else
                error('Unsupported type of compression');
            end
        otherwise
            error('Unsupported option specified');
    end
    
end

