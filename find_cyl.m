function [x0, y0, r, h] = find_cyl(neuropoints)

x0 = mean(neuropoints(:,1));
y0 = mean(neuropoints(:,2));

r = max(eu_dist([x0,y0],neuropoints(:,1:2)));
h = max(neuropoints(:,3));