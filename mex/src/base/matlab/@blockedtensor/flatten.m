function t = flatten(bt)
%FLATTEN flatten a blocked tensor.
%
%   X = FLATTEN(A) creates a tensor ofelements from the blocked tensor A.
%
%   Examples
%   A = blockedtensor({tensor(rand(2,2),[2,2]), tensor(rand(2,2),[2,2])}, [4,2],[2,2])
%   X = flatten(A) %<-- tensor of size 4 x 2 with elements filled
%   appropriately
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

blk_size = bt.block_size;
blked_size = bt.flat_size ./ bt.block_size;
t_size = bt.flat_size;
t_data = zeros(t_size);

nBlks = prod(blked_size);

for i = 1:nBlks
    blkSub = tt_ind2sub(blked_size,i);
    updateregion = arrayfun(@(x,y) {(x-1)*y+1:x*y}, blkSub, blk_size);
    t_data(updateregion{:}) = bt.data_blocks{i}.data;
end

t = tensor(t_data, t_size);
return;