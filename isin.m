function [bool] = isin(x, list)

lenx = length(x);
bool = zeros(lenx,1);

if lenx == 1
    bool = sum(x == list);
else
    for i = 1:lenx
        bool(i) = sum(x(i) == list);
    end
end

bool = logical(bool);