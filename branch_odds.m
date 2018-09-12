function [b_odds] = branch_odds(n_branches)

if n_branches == 5
    b_odds = 0.002; %0.2% of branches are quint    
elseif n_branches == 4
    b_odds = 0.006; %.6% of branches are quad
elseif n_branches == 3
    b_odds = 0.107; %10.7% of branches are tri
else %88.5% of branches are bi
    b_odds = 1; 
end