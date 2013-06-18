function t = test_blockedpsym()
test_blocked_sym_tensor()
test_blocked_psym_tensor_work()
%test_blocked_psym_tensor_break()
return;

function [] = test_blocked_sym_tensor()
%test if symmetric creation works
order = 3;
sym = {1:order};
flat_dimension = 4;
blk_dimension = 2;
blked_dimension = flat_dimension / blk_dimension;

blk_size = genTensorShape(blk_dimension, sym);

uniqueSubscripts = tt_createUniqueIndices(blked_dimension, sym);

index_symmetries = tt_ind2symGroup(uniqueSubscripts);
%Use create_problem to create the blocks with correct symmetry
problems = cell(size(index_symmetries, 1), 1);

for i=1:size(problems)
	problems{i} = create_problem('Size', blk_size, 'Symmetric', block_sym(index_symmetries{i}, sym));
end
data_blks = cellfun(@(x) x.Data, problems, 'UniformOutput', false);
bsymt = blockedpsymtensor(data_blks, order, (flat_dimension), (blk_dimension));

backret = TLA_blockedpsymtensor(bsymt)

return

function [] = test_blocked_psym_tensor_work()
%test if p-symmetric works

order = 5;
psym = {[1 3], [4], [2 5]};
flat_dimension = [4 8 6];
blk_dimension = [4 1 2];
blked_dimension = flat_dimension ./ blk_dimension;

blk_size = genTensorShape(blk_dimension, psym);

uniqueSubscripts = tt_createUniqueIndices(blked_dimension, psym);

index_symmetries = tt_ind2symGroup(uniqueSubscripts);
%Use create_problem to create the blocks with correct symmetry
problems = cell(size(index_symmetries, 1), 1);
for i=1:size(problems)
	problems{i} = create_problem('Size', blk_size, 'Symmetric', block_sym(index_symmetries{i}, psym));
end
data_blks = cellfun(@(x) x.Data, problems, 'UniformOutput', false);
bpsymt = blockedpsymtensor(data_blks, order, flat_dimension, blk_dimension, psym);

%test if correctly rejects bad psym input
disp 'psym_blocks'
for i=1:numel(bpsymt.data_blocks)
    disp(bpsymt.data_blocks{i})
end

psym_backret = TLA_blockedpsymtensor(bpsymt)

return

function [] = test_blocked_psym_tensor_break()
%test if correctly rejects psym input

return

function ten_size = genTensorShape(symGroupDim, psym)
order = numel([psym{:}]);
unordered_size = cell(1,numel(psym));
for i=1:numel(psym)
    unordered_size{i} = repmat(symGroupDim(i), [1 numel(psym{i})]);
end

unordered_size = [unordered_size{:}];
invperm = zeros(1,order);
invperm([psym{:}]) = 1:order;

ten_size = unordered_size(invperm);
return;

function blkSym = block_sym(symIndex, symTensor)
order = numel(symIndex{:});
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


