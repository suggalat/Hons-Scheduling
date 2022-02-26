n = 4;
C = cell(1,n);
[C{:}] = ndgrid([1,-1]);
C = cellfun(@(a)a(:),C,'Uni',0);
M = [C{:}]
