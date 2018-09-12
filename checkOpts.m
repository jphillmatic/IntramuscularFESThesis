function [is_valid] = checkOpts(tx_new,tx_old)

if any(tx_new<0)
    is_valid = 0;
    return
elseif