function [npnew] = subset_tree(neuropoints, nodes)

npnew = neuropoints(ismember(neuropoints(:,6),nodes),:);
