function [MU_nodes] = simple_MU_nodes(nMu, neuropoints)
% Select NMJ nodes which exist only at the highest level of the nerve tree

MU_nodes = {};

ndepth = max(neuropoints(:,4));
nmj_nodes = neuropoints(neuropoints(:,4) == ndepth, 6);
selected_nmj = sort(datasample(nmj_nodes, nMu, 'Replace', false));

% [x0, y0, r, h] = find_cyl(neuropoints);
% std_x = 0;
% std_y = 0;
% inc = 0;
% 
% while (std_x < r/3) | (std_y < r/3)
%     selected_nmj = sort(datasample(nmj_nodes, nMu, 'Replace', false));
%     std_x = std(neuropoints(selected_nmj,1));
%     std_y = std(neuropoints(selected_nmj,2));
%     inc = inc + 1;
%     if inc > 100
%         error('Unable to achieve separation')
%     end
% end

for i = 1:nMu
    p0 = selected_nmj(i);
    MU_nodes{i} = p0;
    t_parent = p0;
    while t_parent > 1
        t_parent = find_parent(MU_nodes{i}(end), neuropoints);
        MU_nodes{i} = horzcat(MU_nodes{i}, t_parent);
    end
end