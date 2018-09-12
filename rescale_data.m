function [rescaledData] = rescale_data(data, window_scales)
rescaledData = [];
for i = 1:size(window_scales,2)
    newData = [];
    for j = 1:size(window_scales,1)
        newData = vertcat(newData, data * window_scales(j));
    end
    rescaledData = horzcat(rescaledData, newData);
end