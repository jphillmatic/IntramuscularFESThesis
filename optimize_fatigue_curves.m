function [max_differences] = optimize_fatigue_curves(mod,fmax1,fatigue,Fsingle,dt,window_len);

%Calculate max in each window for single motor unit model
fmax_single = get_max(Fsingle, window_len/dt);

%Calculate percent decline of max force in each train for single motor unit 
percent_decline = fmax_single ./ fmax_single(1);

%Calculate the percentage of total %decline which happens in each train
p2_decline = (1-percent_decline) / max(1-percent_decline);

%Calculate proportion of of original force of each train for each MU
force_prop = 1 - fatigue.*p2_decline*mod;

%Calculate the new max forces in each window, by taking the product of the
%max force first train of the single MU and force_prop, the proportion of
%single-MU-1st-train-max represetnted by each MU in each train window
f_new = sum(force_prop.*fmax1,2);

%return the mean squared error between the window maxes of single motor
%unit and multi motor unit models
max_differences = d_mse(fmax_single,f_new);
