function [daughters] = find_daughters(nodes, neuropoints)

daughters = neuropoints(logical(sum(neuropoints(:,7) == nodes,2)), 6);

end