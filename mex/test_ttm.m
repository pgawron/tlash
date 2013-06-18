datablks = {tensor([1 2 3 4], [2,1,2]), tensor([2 3 4 5],[2 1 2]), tensor([3 4 5 6],[2 1 2]), tensor([4 5 6 7], [2 1 2])};
bt = blockedtensor(datablks, [4 1 4], [2 1 2]);

datablks = {tensor([1 2 3 4],[2 2]), tensor([1 2 3 4],[2 2]), tensor([1 2 3 4], [2 2]), tensor([1 2 3 4], [2 2])};
bm = blockedtensor(datablks, [4 4], [2 2]);

bC = TLA_ttm(bt, bm, 1);
check = flatten(bC);

A = flatten(bt)
flatB = flatten(bm);
B = reshape(flatB.data, bm.flat_size)
diff = check - ttm(A,B,1);
diff = max(diff.data(:))
