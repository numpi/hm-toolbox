function opt = hmatrixoption(key, value)
%hmatrixOPTION Set or get an option for the hmatrix toolbox.
%
% Valid options are:
%   'block-size': Integer representing the minimum block size.
%   'threshold': Value used for off-diagonal truncation.
%   'compression': String indicating the strategy of compression (qr or svd)

global hmatrix_block_size
global hmatrix_threshold
global hmatrix_compression

if isempty(hmatrix_block_size)
    hmatrix_block_size = 256;
end

if isempty(hmatrix_compression)
    hmatrix_compression = 'svd';
end

if isempty(hmatrix_threshold)
    hmatrix_threshold = 1e-12;
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

