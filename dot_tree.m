function [node_points, dtrs] = dot_tree(neuropoints, nodes)

if nargin == 2
    neuropoints = neuropoints(isin(neuropoints(:,6), nodes),:);
end

% figure(1)
% hold on
parents = unique(neuropoints(:,7));
node_points = [0 0 0];
dtrs = [0];
schwann_len = 10e-4 * 115; %(cm) avg nerve diameter  = 10e-6, scale factor -> node distance  = 115
% schwann_len = 0.1; %cm

for i = 1:length(parents)
    p = parents(i);
    
    pnp = neuropoints(neuropoints(:,6) == p,:);
    
    daughters = find_daughters(p, neuropoints);
    
    
    if ~(isempty(daughters) || isempty(pnp))
        for j = 1:length(daughters)
            d = daughters(j);

            dnp = neuropoints(neuropoints(:,6) == d,:);
            
%             line([pnp(1), dnp(1)], [pnp(2), dnp(2)], [pnp(3), dnp(3)])
            vec_len = pdist2(dnp(1:3), pnp(1:3));
            npoints = vec_len/schwann_len;
            new_pts = linspaceNDim(pnp(1:3), dnp(1:3),max(round(npoints) + 1,3))';
            new_pts = new_pts(2:end,:);
            node_points = vertcat(node_points, new_pts);
            dtrs = vertcat(dtrs, repmat(d,size(new_pts,1),1));
        end
    end
end
dtrs(1) = dtrs(2);
% scatter3(node_points(:,1), node_points(:,2), node_points(:,3), 'o', 'MarkerEdgeAlpha',.2)
