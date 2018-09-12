function [tree_length] = getTreeLength(neuropoints, nodes)

if nargin == 2
    neuropoints = subset_tree(neuropoints, nodes);
end

tree_length = 0;

for i = 1:size(neuropoints,1)
    t_nod = neuropoints(i,6);
    t_dtrs = find_daughters(t_nod, neuropoints);
    
    for j = 1:length(t_dtrs)
        
        t_dist = pdist2(neuropoints(i,1:3), neuropoints(neuropoints(:,6) == t_dtrs(j),1:3));
        
        tree_length = tree_length + t_dist;
    end
end
    

