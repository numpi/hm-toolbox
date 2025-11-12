classdef tel
    %TEL Telescopic decomposition of a (rectangular?) matrix A, represented
    % as A = blkdiag(D) + blkdiag(U) * M * blkdiag(V)', where M is
    % recursively represented as a @tel object. 
    
    properties
        % cell-array of the diagonal blocks at this level, traversed 
        D

        % cell-arrays of left and right generators
        U
        V

        % inner telescopic matrix, defined recursively; at the bottom of
        % the recursion, this is a dense matrix
        M
    end
    
    methods
        function obj = tel(varargin)
            if nargin == 4
                obj.D = varargin{1}; 
                obj.U = varargin{2};
                obj.V = varargin{3};
                obj.M = varargin{4};
            end
        end
    end
end

