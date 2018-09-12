function [t_diffs] = optimize_delay_param(tx, d_ct)
dt = 0.1;
% tx=200;
if length(tx) == 3
    F = calc_MU_force(1, 1, dt, ["t1", "t2","Km_rest"], [tx(1); tx(2); tx(3)],"c100");
elseif length(tx) == 4
    F = calc_MU_force(1, 1, dt, ["t1", "t2","Km_rest", "tc"], [tx(1); tx(2); tx(3);tx(4)],"c100");
else
    F = calc_MU_force(1, 1, dt, ["t1", "t2"], [tx(1); tx(2)],"c100");
end
[mx, indmx] = get_max(F);
indhm = find(F > (mx/2),1,'last');
% indmx*dt
% d_ct
ct_diff = abs(indmx*dt - d_ct);
hm_diff = abs(indhm*dt - (d_ct+20)*2);
t_diffs = ct_diff + hm_diff;
% ct_hrt = [indmx,indhm];
% plot(F)
% hold on