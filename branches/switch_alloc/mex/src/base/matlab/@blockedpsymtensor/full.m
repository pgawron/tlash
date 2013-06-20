function bt = full(t)
%FULL Convert blockedpsymtensor to blockedtensor.
%
%   X = FULL(A) creates a blocked tensor from the blocked psym tensor A. 
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

order = numel(t.flat_size);
bt_flat_size = t.flat_size;
bt_block_size = t.block_size;

bt_blked_size = bt_flat_size ./ bt_block_size;

sym = t.sym;
nBlocks = prod(bt_blked_size);

uniqueIndices = zeros(nBlocks, order);
uniqueIndexPerms = zeros(nBlocks, order);
for i=1:nBlocks
    [uniqueIndexPerm, uniqueIndex] = getUniqueIndex(tt_ind2sub(bt_blked_size, i), sym);
    uniqueIndexPerms(i,:) = uniqueIndexPerm;
    uniqueIndices(i,:) = uniqueIndex;
end
%uniqueIndices = uniqueIndices;
%uniqueIndexPerms = uniqueIndexPerms;
%merged = sortrows([uniqueIndices, uniqueIndexPerms], order:-1:1);
%uniqueIndices = merged(1:order);
%uniqueIndexPerms = merged(order+1:end);
symIndices = sortrows(unique(uniqueIndices, 'rows'), order:-1:1);

data_blks = cell(1,nBlocks);
for i = 1:nBlocks
    uniqueIndex = uniqueIndices(i,:);
    [~, uniqueLinIndex] = ismember(uniqueIndex, symIndices, 'rows');
    data_blks{i} = ipermute(t.data_blocks{uniqueLinIndex}, uniqueIndexPerms(i,:));
end
bt = blockedtensor(data_blks, bt_flat_size, bt_block_size);

%disp 'psymtensor'
%for i = 1:size(symIndices,1)
%    disp(t.data_blocks{i})
%end
%disp 'blockedtensor'
%for i = 1:nBlocks
%    disp(bt.data_blocks{i})
%end
end

function [perm, ind] = getUniqueIndex(blk_index, sym)
    order = numel(blk_index);
    ind = blk_index;
    perm = 1:order;
    for i = 1:numel(sym)
        symGroup = sym{i};
        sorted = sortrows([ind(symGroup); symGroup]', 1);
        ind(symGroup) = (sorted(:,1))';
        perm(symGroup) = (sorted(:,2))';
    end
end