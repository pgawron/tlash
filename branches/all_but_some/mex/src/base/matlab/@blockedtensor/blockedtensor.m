function t = blockedtensor(varargin)
%Remove FLATSIZE/BLKSIZE redundancy.  Can get away with just the size of
%the tensor?
%BLOCKEDTENSOR Create blocked tensor.
%
%   X = BLOCKEDTENSOR({A},FLATSIZE, BLKSIZE) creates a blocked tensor from the multidimensional
%   array of tensors A. The FLATSIZE and BLKSIZE arguments specify the desired shape of A.
%
%   X = BLOCKEDTENSOR(S) copies a blockedtensor S.
%
%   Examples
%   X = blockedtensor({tensor([...],[2,1,4]), tensor([...],[2,1,4]), ...}, [2,4,8],[2,1,4]) %<-- Blocked Tensor of
%   size 6 x 4 x 8 with block size 2 x 1 x 4
%
%   See also TENSOR/NDIMS.
%
%MATLAB Tensor Toolbox.
%Copyright 2012, Sandia Corporation.

% This is the MATLAB Tensor Toolbox by T. Kolda, B. Bader, and others.
% http://www.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2012) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in the file LICENSE.txt


% EMPTY/DEFAULT CONSTRUCTOR
if nargin == 0
	t.data_blocks = {};
    t.flat_size = [];
    t.block_size = [];
    t = class(t, 'blockedtensor');
    return;
end

% CONVERSION/COPY CONSTRUCTORS
% Note that we pass through this if/switch statement if the first argument
% is not any of these cases.
if (nargin == 1)
    v = varargin{1};
    switch class(v)
        case 'blockedtensor',
            %COPY CONSTRUCTOR
            t.data_blocks = v.data_blocks;
            t.flat_size = v.flat_size;
            t.block_size = v.block_size;
            t = class(t, 'blockedtensor');
            return;
        case 'blockedpsymtensor'
            t = full(v);
            return;
%        case 'tensor',   
%            t.data = v.data;
%            t.size = v.size;
%            t = class(t, 'tensor');
%            return;
    end
end

% CONVERT A MULTIDIMENSIONAL ARRAY
if (nargin == 3)
    flat_size = varargin{2};
    if ndims(flat_size) ~= 2 || size(flat_size,1) ~= 1
            error('Second argument must be a row vector.');
    end
    
    block_size = varargin{3};
    if ndims(block_size) ~= 2 || size(block_size,1) ~= 1
            error('Third argument must be a row vector.');
    end
    
    if(size(flat_size, 2) ~= size(block_size, 2))
        error('flat_size and blk_size must be same order');
    end
    
    blked_size = flat_size ./ block_size;
    order = numel(flat_size);
    numBlocks = prod(blked_size);
    
    dataBlks = varargin{1};
    if ndims(dataBlks) ~= 2 || size(dataBlks,1) ~= 1
            error('First argument must be a cell');
    end 
    
    if numel(dataBlks) ~= numBlocks
        error('First argument has incorrect number of blocks')
    end
    
    if(any(cellfun(@(x) size(x.size,2), dataBlks) ~= order))
        error('All data blocks must be same order tensor as "order"')
    end
    
    dataSizes = cellfun(@size, dataBlks, 'UniformOutput', false);
    guard = cellfun(@(x) any(x ~= block_size), dataSizes);
    if(any(guard))
        error('All data blocks must have size equal to blk_size')
    end
    
    %create the blocked tensor
    t.flat_size = flat_size;
    t.block_size = block_size;
    t.data_blocks = dataBlks;
    t = class(t, 'blockedtensor');
    
    return;
end


error('Unsupported use of function TENSOR.');