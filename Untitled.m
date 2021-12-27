load bnodes
[edges_cell, nodes_cell, seeds] = voronoi_bbx(bnodes, 256, 128);

figure; scatter(seeds(:,1),seeds(:,2)); hold on;
figure
for ele = 1:size(nodes_cell)
    nodes = nodes_cell{ele};
    scatter(nodes(:,1),nodes(:,2)); hold on;
end