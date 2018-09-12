function [npnew] = removeOutliers(neuropoints)

[x0, y0, r, h] = find_cyl(neuropoints);
p_outside = neuropoints(sqrt(sum(neuropoints(:,1:2).^2,2)) > r*0.75,6);

for i = 1:length(neuropoints)
    tpp = parent_path(neuropoints(i,6), neuropoints);
    if any(ismember(tpp,p_outside))
        p_outside(end+1) = neuropoints(i,6);
    end
end


npnew = neuropoints(~ismember(neuropoints(:,6),p_outside),:);

nodes_new = 1:size(npnew,1);

m = containers.Map(npnew(:,6),nodes_new);

npnew(:,6) = nodes_new;
for i = 2:size(npnew,1)
    npnew(i,7) = m(npnew(i,7));
end