function [tA, mB, tC] = test_sttsm(order, flatSizeA, flatSizeC, blkSizeA, blkSizeC)
	A = create_random_blocked_sym_tensor(order, flatSizeA * ones(1,order), blkSizeA * ones(1,order));
	B = create_random_blocked_tensor([flatSizeC, flatSizeA], [blkSizeC, blkSizeA]);	
	C = TLA_sttsm(A, B);
    disp(C)
    %Check answer
    
    checkA = tensor(A);
    
    fb = flatten(B);
    checkB = fb.data;
    disp 'checking C'
    for i = 1:numel(C.data_blocks)
	fprintf('matlab blk %d: ', i);
	C.data_blocks{i}
    end
    checkC = tensor(C);
    
    tA = checkA;
    mB = checkB;
    tC = checkC;
    bmult = cell(1, order);
    bmult(:) = {checkB};
    checkRes = ttm(checkA, bmult);
    diff = norm(checkC - checkRes)
return;
end

function bpsymt = create_blocked_psym_tensor(order, flat_size, block_size, sym)

flat_symGroup_size = cellfun(@(x) flat_size(x(1)), sym);
block_symGroup_size = cellfun(@(x) block_size(x(1)), sym);
blocked_symGroup_size = flat_symGroup_size ./ block_symGroup_size;

uniqueSubscripts = tt_createUniqueIndices(blocked_symGroup_size, sym);

index_symmetries = tt_ind2symGroup(uniqueSubscripts);
%Use create_problem to create the blocks with correct symmetry
data_blks = cell(size(index_symmetries, 1), 1);

for i=1:size(data_blks)
	data_blks{i} = create_block(block_size, block_sym(index_symmetries{i}, sym));
end

bpsymt = blockedpsymtensor(data_blks, order, flat_symGroup_size, block_symGroup_size);

return
end

function bsymt = create_random_blocked_sym_tensor(order, flat_size, block_size)
sym = {1:order};
bsymt = create_blocked_psym_tensor(order, flat_size, block_size, sym);
return
end

function bt = create_random_blocked_tensor(flat_size, block_size)
nBlocks = prod(flat_size ./block_size);
data_blocks = cell(1, nBlocks);

for i=1:numel(data_blocks)
    data_blocks{i} = tensor(rand(block_size));
end

bt = blockedtensor(data_blocks, flat_size, block_size);
return;
end

function blk = create_block(block_size, sym)
blk = tensor(1, ones(1, numel(block_size)));

for i=1:numel(sym)
    symGroup = sym{i};
    vec = rand(block_size(i), 1);
    bmult = cell(1,numel(symGroup));
    bmult(:) = {vec};
    blk = ttm(blk, bmult, symGroup);
end
end

function blkSym = block_sym(symIndex, symTensor)
order = numel([symIndex{:}]);
toFind = 1:order;

blkSym = {};
while ~isempty(toFind)
		modeToFind = toFind(1);
		
		indexSymGroup = symIndex{cellfun(@(x) any(x == modeToFind), symIndex)};
		tensorSymGroup = symTensor{cellfun(@(x) any(x == modeToFind), symTensor)};

		thisSym = indexSymGroup(ismember(indexSymGroup, tensorSymGroup));
		toFind = toFind(~ismember(toFind, thisSym));

		blkSym = cat(2,blkSym, thisSym);
end
return;
end
