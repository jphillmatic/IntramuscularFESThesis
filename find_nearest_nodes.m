function [npnew] = find_nearest_nodes(nmj_pos, neuropoints)
nnodes = zeros(size(nmj_pos,1),4);
numnodes = max(neuropoints(:,6));
[a, ind] = pdist2(neuropoints(:,1:3),nmj_pos,'euclidean','Smallest',1);

for i = 1:length(ind)
    nnodes(i,1) = neuropoints(ind(i),4)+1;
    nnodes(i,4) = neuropoints(ind(i),6);
    nnodes(i,3) = numnodes + i;
    nnodes(i,2) = length(find_daughters(neuropoints(ind(i),6),neuropoints))+1;
end
    
nnodes = horzcat(nmj_pos, nnodes);

npnew = vertcat(neuropoints, nnodes);

pnodes = [];

for i = 1:size(nmj_pos,1)
    pnodes  = vertcat(pnodes,parent_path(nnodes(i,6),npnew));
end

npnew = subset_tree(npnew, unique(pnodes));
