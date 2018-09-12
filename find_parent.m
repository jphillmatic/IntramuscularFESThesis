function [parent] = find_parent(node , neuropoints)

parent = neuropoints((neuropoints(:,6) == node),7);

end