function [nmj_nodes] = find_nmj(neuropoints)
max_dist = max(neuropoints(:,4));

nmj_nodes =[];
for i = 1:size(neuropoints,1)
    ndn = find_daughters(neuropoints(i,6), neuropoints);
    if isempty(ndn) %&& (neuropoints(i,4) > max_dist/4)
        nmj_nodes = vertcat(nmj_nodes, neuropoints(i,6));
    end
end
    