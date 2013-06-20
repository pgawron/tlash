function symGroup = tt_ind2symGroup(ind)

%IND2SYMGROUP Create the symmetries associated with supplied index
%
%   SYM = TT_IND2SYMGROUP(IND) creates a symGroup reprsenting the symmetry associated with the index IN
%
%   Examples
%   SYM = tt_ind2symGroup([1 1 2 3 2 4]) %<-- Symgroup {[1 2],[3 5],[4],[6]}
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

uniqueVals = unique(ind);
symGroup = arrayfun(@(x) find(ind == x), uniqueVals, 'UniformOutput', false);
return;
