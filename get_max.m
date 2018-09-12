function [maxval, ind] = get_max(data, window_size)

if nargin == 1
    window_size = length(data);
end

n = size(data,2);
n_windows = floor(length(data) / window_size);
maxval = zeros(n_windows,n);
ind = zeros(n_windows,n);

for mu = 1:n
    for i = 1:n_windows
        frange = window_size*(i-1)+1:window_size*(i);
        [maxval(i,mu), ind(i,mu)] = max(data(frange,mu));
        ind(i,mu) = ind(i,mu) + window_size * (i-1);
    end
end