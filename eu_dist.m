function [distances] = eu_dist(p_electrode, p_nodes)

p_nodes = p_nodes - p_electrode;

distances = sqrt(sum(p_nodes.^2,2));