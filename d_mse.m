function [ms_e] = d_mse(seq1, seq2)

ms_e = sum(sqrt((seq1 - seq2).^2))/length(seq1);