function [msee] = optimize_fatigue_param(mfat, P, ForceMaxes, dt, trains, window_len)

dt = 1; %ms

trains = 75;
window_len = 1750;

Force = calc_MU_force(1, trains, dt, ["A_rest", "t_fat"], [P; mfat]);

newMax = get_max(Force, window_len/dt);

msee = d_mse(ForceMaxes, newMax);
