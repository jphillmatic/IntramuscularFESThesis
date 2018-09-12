function [m_ct] = optimize_delay_time(mct, P, ForceMaxes, dt, trains, window_len)
%Optimize delay time so that the peak force is coincident with Ding model
%Time to peak (P) and max twitch force(T) from Fugelvand 1993  
forceRange = 120; % 10-fold force diffecence between largest and smallest MU
fatigueRange = 180; % range of fatigue rates across the motor units (300 best)
contractionRange = 3; % Range in contraction times,
ct_max = 90;
% FatFac = 0.0225;    % fatigue factor (FF/S) percent of peak force of MU per second 
tick_breaks = 3;
[P, fatigue, arborLengths, m_csa, ct]  = define_motor_pool(nMu, forceRange, fatigueRange, contractionRange);

% param_names = ["A_rest" "t_fat" "t1" "t2" "aA" "aT1" "aKm"];
% param_vals = [5.1 53.4*1000 43.8 124.4 -8.8e-7 5.7e-6 1.9e-8];
F = calc_MU_force(nMu, 1, dt, ["t1" "t2"], [ct; ct/0.3],'c100');
F0 = calc_MU_force(1, 1, dt, [], [], 'c100')
figure(3)
hold on
clf
[bmax, ibmax] = get_max(F);
for pmu = 1:nMu
    mu = nMu +1 - pmu;
    % sum(force_scales,2);
    plot(F(:,mu), 'Color', scale_gradient(pmu,:))
    hold on
    plot(ibmax(mu),bmax(mu),'ro')
end
plot(F0, 'Color', reed)
plot(sum(F,2), 'Color', bloo)

