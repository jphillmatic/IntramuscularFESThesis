function [MU_nodes] = select_MU_nodes(arborLengths, neuropoints)

%Cell structure, each cell is unique motor unit, and each contains all
%nodes of the motor unit

MU_nodes = {};
nMu = length(arborLengths);

% ndepth = max(neuropoints(:,4));
% nmj_points = neuropoints(neuropoints(:,4) == ndepth, 6);
% nmj_points = datasample(nmj_points, nMu, 'Replace', false);
    
nmj_points = find_nmj(neuropoints);
nmj_sample = datasample(nmj_points, nMu, 'Replace', false);

%partitioned muscle
partitioned = true;
if partitioned
    branch_start_nodes = neuropoints(neuropoints(:,4)==2,6);
    rpmt = repmat([1:length(branch_start_nodes)],ceil(nMu/length(branch_start_nodes)),1);
    branch_groups = datasample(rpmt(:), nMu, 'Replace', false);
    
    for i = 1:length(nmj_points)
        pp(i) = find(ismember(branch_start_nodes, parent_path(nmj_points(i),neuropoints)));
    end
    
end

for i = 1:nMu
    p0 = nmj_sample(i);
    if partitioned
        p0 = datasample(nmj_points(pp == branch_groups(i)),1);
    end
    MU_nodes{i} = p0;
    t_parent = p0;
    while t_parent > 1
        t_parent = find_parent(MU_nodes{i}(end), neuropoints);
        MU_nodes{i} = horzcat(MU_nodes{i}, t_parent);
    end
    
    %temporary arbor length
%     al = length(MU_nodes{i});
    al = getTreeLength(neuropoints, MU_nodes{i});
    
    while al < arborLengths(i)
% ADD NODES BY NMJ        
        nnmj = datasample(nmj_points,1);
        if partitioned
            nnmj = datasample(nmj_points(pp == branch_groups(i)),1);
        end
        nd = nnmj;
        t_parent = nnmj;
        while t_parent > 1
            t_parent = find_parent(nd(end), neuropoints);
            nd = horzcat(nd, t_parent);
        end

% % ADD NODES BY INDIVIDUAL DAUGHTER
%         t_daughters = unique(find_daughters(MU_nodes{i}, neuropoints));
%         t_daughters = t_daughters(~isin(t_daughters, MU_nodes{i}));
%         
%         %select a daughter
%         %nd = datasample(t_daughters,1);
%         valid_daughter = 0;
%         while ~valid_daughter
%             nd = datasample(t_daughters,1);
%             p_temp = find_parent(nd,neuropoints);
%             n_competing_daughters = sum(ismember(find_daughters(p_temp,neuropoints), MU_nodes{i}));
%             
%             branch_rand = rand(1);
%             
%             if branch_rand <= branch_odds(n_competing_daughters)
%                 valid_daughter = 1;
%             end            
%         end
        
        %add new node to the motor unit
        MU_nodes{i} = unique(horzcat(MU_nodes{i}, nd));
        %recalculate arbor length
%         al = length(MU_nodes{i});
        al = getTreeLength(neuropoints, MU_nodes{i});
    end
end