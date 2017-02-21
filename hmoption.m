function opt = hmoption(key, value)
%HMOPTION Set or get an option for the hm toolbox. 
%
% Valid options are: 
%   'block-size': Integer representing the minimum block size. 

global hm_block_size

if isempty(hm_block_size)
	hm_block_size = 96;
end

if ~exist('value', 'var')
	switch key
		case 'block-size'
			opt = hm_block_size;
		otherwise
			error('Unsupported option specified');
	end
else
	switch key
		case 'block-size'
			hm_block_size = value;
		otherwise
			error('Unsupported option specified');
	end

end

