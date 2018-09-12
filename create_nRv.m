function [node_points] = create_nRv(neuropoints, nodes)

if nargin == 2
    neuropoints = neuropoints(isin(neuropoints(:,6), nodes),:);
end


parents = unique(neuropoints(:,7));
% node_points = zeros(100 * length(neuropoints(:,6)),4);
node_points = [];
schwann_len = 0.115; %cm

for i = 1:length(parents)
    p = parents(i);
    
    pnp = neuropoints(neuropoints(:,6) == p,:);
    
    daughters = find_daughters(p, neuropoints);
    
    new_pts = [];
    
    if ~(isempty(daughters) || isempty(pnp))
        for j = 1:length(daughters)
            d = daughters(j);

            dnp = neuropoints(neuropoints(:,6) == d,:);
            
%             line([pnp(1), dnp(1)], [pnp(2), dnp(2)], [pnp(3), dnp(3)])
            vec = dnp(1:3) - pnp(1:3);
            vec_len = sqrt(sum(vec.^2));
            npoints = vec_len/schwann_len;
            new_pts = linspaceNDim(pnp(1:3), dnp(1:3),round(npoints) + 1)';
            
            node_points = vertcat(node_points,[new_pts,repmat(d,length(new_pts),1)]);
        end
    end
end
% figure(1)
% hold on
% scatter3(node_points(:,1), node_points(:,2), node_points(:,3), '.', 'MarkerEdgeAlpha',.2)
