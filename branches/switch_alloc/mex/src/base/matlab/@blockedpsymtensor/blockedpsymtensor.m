function t = blockedpsymtensor(varargin)
%Remove FLATSIZE/BLKSIZE redundancy.  Can get away with just the size of
%the tensor?
%BLOCKEDPSYMTENSOR Create blocked (partially)symmetric tensor.
%
%   X = BLOCKEDPSYMTENSOR({A}, ORDER, FLATSIZE, BLKSIZE) creates a blocked, 
%       symmetric, order-ORDER tensor from the multidimensional
%       array of tensors A. The FLATSIZE and BLKSIZE argument specifies 
%       the desired shape of symmetric groups of A.
%
%   X = BLOCKEDPSYMTENSOR({A}, ORDER, FLATSIZE, BLKSIZE, SYM) creates a blocked, 
%       SYM-symmetric, order-ORDER tensor from the multidimensional
%       array of tensors A. The SIZE argument specifies the desired 
%       shape of symmetric groups of A.
%
%   X = BLOCKEDTENSOR(S) copies a blockedtensor S.
%
%   Examples
%   SYM = {[1 2 3],[4 5]}
%   X = blockedpsymtensor({tensor(...),...}, 5, [4 6], [2 3], SYM) %<-- an
%   order-5 ((1 2 3) (4 5))-symmetric 4 x 4 x 4 x 6 x 6 tensor 
%   with 2 x 2 x 2 x 3 x 3 blocks
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
    t.sym = {};
    t = class(t, 'blockedpsymtensor');
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
        case 'blockedpsymtensor',
            t.data_blocks = v.data_blocks;
            t.flat_size = v.flat_size;
            t.block_size = v.block_size;
            t.sym = v.sym;
            t = class(t, 'blockedpsymtensor');
    end
end

% CONVERT A MULTIDIMENSIONAL ARRAY
if (nargin == 4)
    t = blockedpsymtensor(varargin{:}, {1:varargin{2}});
    return;
elseif (nargin == 5)
    data_blks = varargin{1};
    order = varargin{2};
    symGroupSizes = varargin{3};
    symGroupBlkSizes = varargin{4};
    sym = varargin{5};
    
    if ndims(symGroupSizes) ~= 2 || size(symGroupSizes,1) ~= 1
        error('Third argument must be a row vector.');
    end
    
    if ndims(symGroupBlkSizes) ~= 2 || size(symGroupBlkSizes,1) ~= 1
        error('Fourth argument must be a row vector.');
    end
    
    if any(mod(symGroupSizes, symGroupBlkSizes) ~= 0)
        error('blkSize must evenly divide symGroupSizes')
    end
    
    symGroupBlkSizes
    sym
    if(numel(symGroupBlkSizes) ~= numel(sym))
        error('blk_size must have same number of elements as symmetric groups')
    end
    
    if(size(symGroupSizes, 2) ~= size(symGroupBlkSizes, 2))
        error('flat_size and blk_size must be same order');
    end
    
    symGroupBlkedSizes = symGroupSizes ./ symGroupBlkSizes;
    
    flat_size = zeros(1,order);
    block_size = zeros(1,order);
    symGroupLens = cellfun(@(x) size(x,2), sym);
    
    if sum(symGroupLens) ~= order
        error('Sym must involve "order" elements')
    end
    
    symModes = [sym{:}];
    
    if numel(unique(symModes)) ~= order
        error('Sym must contain unique entries')
    end
    
    if numel(unique(symModes)) ~= order || any(symModes > order) || any(symModes < 1)
        error('Sym can only have values 1 <= x <= order')
    end
    
    %think about this a bit harder
    for i=1:numel(sym)
        flat_size(sym{i}) = symGroupSizes(i);
        block_size(sym{i}) = symGroupBlkSizes(i);
    end
    
    nSymGroupUniqueBlocks = cellfun(@(x, y) nchoosek(numel(x) + y - 1, numel(x)), sym, num2cell(symGroupBlkedSizes));
    nUniqueBlocks = prod(nSymGroupUniqueBlocks);
    
    if numel(data_blks) ~= nUniqueBlocks
        error('Incorrect number of data blocks')
    end
    
    
    dataSizes = cellfun(@size, data_blks, 'UniformOutput', false);
    guard = cellfun(@(x) any(x ~= block_size), dataSizes);
    if(any(guard))
        error('All data blocks must have size equal to blk_size')
    end
    
    %check for all data blocks being tensors
    
    %create the blocked psym tensor
    t.flat_size = flat_size;
    t.block_size = block_size;
    t.data_blocks = data_blks;
    t.sym = sym;
    t = class(t, 'blockedpsymtensor');
    
    return;
end


error('Unsupported use of function TENSOR.');