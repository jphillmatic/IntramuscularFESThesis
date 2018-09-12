function [adj_mat] = create_adj_matrix(dtrs, neuropoints)

[udtrs, start_nodes, t2] = unique(dtrs,'first');
[t, last_nodes, t2] = unique(dtrs,'last');

adj_mat = zeros(length(dtrs));
adj_mat(1,1:2) = [-1 1];

for i = 2:length(dtrs)
    %find point towards base
    if ismember(i,start_nodes)
        pt = find_parent(dtrs(i), neuropoints);
        prev_pt = last_nodes(pt == udtrs);         
    else
        prev_pt = i-1;
    end
    adj_mat(i,prev_pt) = 1;
    
    %find points toward nmj
    if ismember(i,last_nodes)
        [t, dtri] = ismember(find_daughters(dtrs(i), neuropoints),dtrs); 
        if any(dtri)
            dt = dtrs(dtri(t));
            for j = 1:length(dt)
                next_pt = start_nodes(dt(j) == udtrs);
                adj_mat(i,next_pt) = 1;
            end
        end
    else
        next_pt = i+1;
        adj_mat(i,next_pt) = 1;
    end
    
    %set current point as negative sum of other connections
    adj_mat(i,i) = -sum(adj_mat(i,:));
    
end
    
    