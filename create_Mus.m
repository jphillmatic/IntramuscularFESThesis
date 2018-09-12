function [MU_nodes] = select_MU_nodes(arborLengths, neuropoints)

%Cell structure, each cell is unique motor unit, and each 
MU_nodes = {};

nMu = length(arborLengths);
ndepth = max(neuropoints(:,4));
top_points = neuropoints(neuropoints(:,4) == ndepth, 6);
top_points = datasample(top_points, nMu);

for i = 1:nMu
    p0 = top_points(i);
    MU_nodes{i} = p0;
    t_parent = p0;
    while t_parent > 1
        t_parent = find_parent(MU_nodes{i}(end), neuropoints);
        MU_nodes{i} = horzcat(MU_nodes{i}, t_parent);
    end
    
    %temporary arbor length
    al = length(MU_nodes{i});
    
    while al < arborLengths(i)
        %find daughters
        t_daughters = unique(find_daughters(MU_nodes{i}, neuropoints));
        t_daughters = t_daughters(~isin(t_daughters, MU_nodes{i}));
        %select a daughter
        nd = datasample(t_daughters,1);
        %add it to he motor unit
        MU_nodes{i} = horzcat(MU_nodes{i}, nd);
        %recalculate arbor length
        al = length(MU_nodes{i});
        
    end
end