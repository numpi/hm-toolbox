function opt = hmoption(key, value)
%HMOPTION Set or get an option for the hm toolbox. 
%
% Valid options are: 
%   'block-size': Integer representing the minimum block size. 
%   'threshold': Value used for off-diagonal truncation. 

global hm_block_size
global hm_threshold

if isempty(hm_block_size)
	hm_block_size = 256;
end

if isempty(hm_threshold)
	hm_threshold = 1e-12;
end

if ~exist('key', 'var')
	error('Please specify a key');
end

if ~exist('value', 'var')
	switch key
		case 'block-size'
			opt = hm_block_size;
		case 'threshold'
			opt = hm_threshold;
		otherwise
			error('Unsupported option specified');
	end
else
	switch key
		case 'block-size'
			if value <= 2
				error('minimum block size must be at least 3');
			else
				hm_block_size = value;
			end
		case 'threshold'
			if value < 0
				error('threshold has to be positive');
			else
				hm_threshold = max(eps, value);
			end
		otherwise
			error('Unsupported option specified');
	end

end

