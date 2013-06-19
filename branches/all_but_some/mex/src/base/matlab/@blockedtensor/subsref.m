function a = subsref(t,s)
%SUBSREF Subscripted reference for a symtensor.
%
%   Examples
%   X = symtensor(1:20, [3,4]);
%   X.data returns the data array ([1:20]).
%   X.size returns the size of the tensor [4,4,4].
%   X(2,3,1) returns that single element of A.
%   X(2,3,:) as of now undefined (TODO)
%
%   See also SYMTENSOR.
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


switch s(1).type    
    case '.'
        switch s(1).subs
            case 'data_blocks'
                a = tt_subsubsref(t.data_blocks,s);
            case 'block_size'
                a = tt_subsubsref(t.block_size,s);
            case 'flat_size'
                a = tt_subsubsref(t.flat_size,s);
            otherwise
                error(['No such field: ', s(1).subs]);
        end
    case '()'
        error('No such functionality () for blockedtensor');
    case '{}'
        error('No such functionality {} for blockedtensor');
    otherwise
        error('Invalid subsref.');
end
