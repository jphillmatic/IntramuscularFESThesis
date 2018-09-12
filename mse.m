function [ms_e] = d_mse(seq1, seq2)

if length(seq1) ~= length(seq2)
    error('Error. Data must have the same length.')
end

ms_e = sum((seq1 - seq1).^2);