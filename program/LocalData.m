cell = [4 1 4 1 3 3 4 1];
cellEnd = zeros(1,5)
NewList = zeros(1,8)
for i = 1:8
    Icell = cell(i);
    cellEnd(Icell+1:end) = cellEnd(Icell+1:end) + 1;
    
end

for i = 8:-1:1
    Icell = cell(i);
    Newlist(cellEnd(Icell+1)) = i;
    cellEnd(Icell+1) = cellEnd(Icell+1) - 1;
end

Newlist
cellEnd