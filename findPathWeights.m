function [mu_thru] = findPathWeights(neuropoints, mu_tree)

nodes = unique(neuropoints(:,6));
mu_thru = zeros(length(nodes),2);

for i = 1:length(nodes)
    t_node = nodes(i);
    mu_thru(i,1) = t_node;
    
    n_pas = 0;
    for j = 1:length(mu_tree)
        if ismember(t_node,mu_tree{j})
            n_pas = n_pas + 1;
        end
    end
    mu_thru(i,2) = n_pas;
end
        