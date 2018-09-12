function [main_branches, n_main_branches] = main_branch_nmjs(neuropoints)

nmj_nodes  = find_nmj(neuropoints);
main_branches = neuropoints(neuropoints(:,4) == 2,6);
n_main_branches = zeros(length(main_branches),1);

for i = 1:length(nmj_nodes)
    pp = parent_path(nmj_nodes(i), neuropoints);
    n_main_branches = n_main_branches + ismember(main_branches,pp);
end
    