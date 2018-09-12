function [parents] = parent_path(node , neuropoints)

parents = zeros(neuropoints(neuropoints(:,6) == node,4),1);
parents(1) = node;
tparent = node;
i = 2;
while tparent > 0
    tparent = find_parent(tparent, neuropoints);
    parents(i) = tparent;
    i = i+1;
end