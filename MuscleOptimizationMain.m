%% Optimization Fatigue Curves 
%Initialize Motor Unit and Stimulation Properties

nMu = 10; %number motor units

LoadColorScales;
global scale_gradient
global scale_discrete
global bloo
global reed

window_len = 1750;%length of a stimulation train (ms)
dt = 0.1; %time step (ms)
forceRange = 120; % 120-fold force diffecence between largest and smallest MU
fatigueRange = 180; % range of fatigue rates across the motor units (300 best)
contractionRange = 2; % range of time to peack twitch force

%Define normalized motor unit properties based on above ranges
[P_twitch, fatigue, arborLengths, m_csa, ct]  = define_motor_pool(nMu, forceRange, fatigueRange, contractionRange);
P_twitch = P_twitch / sum(P_twitch) * 5.1;

%Calculate force output using Ding parameters
[F_Ding03,time] = calc_MU_force(1, 75, dt);
[fmaxD03, imaxD03] = get_max(F_Ding03, 1750/dt);

% %% Parameterize contraction time and rescale P
% Fc100 = calc_MU_force(10, 1, dt, [], [],"c100");
% tx_new = [];
% [mc100 imc100] = get_max(Fc100);
% tx0 = [43.8,124.4,0.3, 20];
% ct_hrt = zeros(nMu,2);
% options = optimoptions('fmincon','Algorithm','interior-point');
% for mu = 1:nMu
%     d_ct = ct(mu);
%     fun = @(tx)optimize_delay_param(tx, d_ct);
% %     valid_minima = false;
% %     tries = 0;
% %     while valid_minima&(tries<10)
%         tx_new(mu,:) = fminsearch(fun, tx0);
% %         valid_minima = checkOpts(tx_new(mu,:), tx_new(mu-1,:));
% %         tries = tries +1;
% %     end
% 
% %     tx_new(mu,:) = fmincon(fun, tx0, [], [], [], [], [0, 0, 0, 0], [],[],options);
% 
% %     optv(mu) = optimize_delay_param(tx_new(mu), d_ct);
% %     Ft = calc_MU_force(1, 1, dt, ["t1", "t2"], [tx_new(mu); tx_new(mu)*3],"c100");
% %     mxt(mu) = get_max(Ft);
% end
% tx_new
% ct_hrt/1000
% %%
% % P_scaling = P_twitch .* mxc100 ./ mxt
% % tx_new
% % optv
% srfc =[];
% for i = 1:8
%     for j= 1:8
%         for k=1:8
%         scale = [i j k]/10;
%             srfc(i,j,k) = optimize_delay_param(tx0.*scale, d_ct);
%         end
%     end
% end
% figure(11)
% clf
% surf(srfc)
% %%
% [Fopt,time] = calc_MU_force(nMu, 1, dt, ["t1", "t2", "Km_rest", "tc"], [tx_new(:,1); tx_new(:,2); tx_new(:,3); tx_new(:,4)],"c100");
% % [Fopt,time] = calc_MU_force(1, 1, dt, ["t1", "t2", "Km_rest"], [tx0(:,1); tx0(:,2); tx0(:,3)],"c100");
% % [Fopt,time] = calc_MU_force(1, 1, dt, ["tc"], 40,"c100");
% 
% [mc100; imc100]
% figure(10)
% clf
% hold on
% [omax, oind] = get_max(Fopt)
% for mu = 1:nMu
%     plot(time/1000, Fopt(:,mu), 'color', scale_gradient(mu,:))
%     
% end
% plot(time/1000,sum(Fc100,2),'-','color',bloo)
% plot(time/1000,sum(Fopt,2),'--','color',reed)

%% Calculate linear scaling factor for fatigue values
%Calculate initial forces generated during 1 train for nMu
ForceT1 = calc_MU_force(nMu, 1, dt);

fmax = zeros(1,nMu);
indmax = zeros(1,nMu);
for mu = 1:nMu
[fmax(:,mu), indmax(:,mu)] = get_max(ForceT1(:,mu), window_len/dt);
end

%Initialize the guess for fatigue scaling factor

%% Optimize fatigue param to match total decline
mod0 = 0.9*nMu;
fun = @(mod)optimize_fatigue_curves(mod,fmax(1,:),fatigue,F_Ding03,dt,1750);
mod = fminsearch(fun, mod0)
%% Generate force scalings
%Calculate max in each window

[fmax_single, indmax_single] = get_max(F_Ding03, window_len/dt);
percent_decline = fmax_single ./ fmax_single(1);
p2_decline = (1-percent_decline) / max(1-percent_decline);

force_prop = 1 - fatigue.*p2_decline*mod;
fm0 = fmax(1,:);
fmax_multi = sum(force_prop.*fm0,2);
mu_pd = force_prop(end,:)./force_prop(1,:);

percent_drop = 1-mu_pd;

d_mse(fmax_single, fmax_multi)
figure(6)
clf
hold on
plot(fmax_single, 'b--')
plot(fmax_multi,'ko')
axis([0 inf 0 inf])
%%
global scale_gradient
global scale_discrete
figure(7)
clf
hold on
F_rescaled = [];

for pmu = 1:nMu
    mu = nMu +1 - pmu;
    F_rescaled(:,mu) = rescale_data(ForceT1(:,mu), force_prop(:,mu));
    % sum(force_scales,2);
    plot(time/1000,F_rescaled(:,mu), 'Color', scale_gradient(pmu,:))
end
plot(time/1000,sum(F_rescaled,2),'Color',scale_discrete(2,:))
size(F_rescaled)
%% Parameterize Fatigue values
%Time Consuming ~2 minutes per motor unit
tic
dt_fatopt = 1;
trains_fatopt = 75;
fmax_mu = fm0.*force_prop;
mfat = zeros(1,nMu);
mfat0 = 1517;
for mu = 1:nMu
    fun = @(ofat)optimize_fatigue_param(ofat, P_twitch(mu), fmax_mu(:,mu), dt_fatopt, trains_fatopt, window_len);
    mfat(mu) = fminsearch(fun, mfat0)
%         mfat(mu) = fmincon(fun, mfat0, [], [])

%     mfat0 = mfat(mu);
end
toc
%% Calculate muscle force using optimized values
muscle_params = [P_twitch; mfat];
[Fnew, time] = calc_MU_force(nMu, 75, 0.1, ["A_rest", "t_fat"], muscle_params);
%% Plot motor unit curves
figure(8)
clf
hold on
plot(F_Ding03,'--','Color',scale_discrete(1,:))
% plot(sum(Fnew,2),'--','Color',scale_discrete(nMu+3,:))
plot(sum(Fnew,2),'Color',scale_discrete(2,:))
for pmu = 1:nMu
    mu = nMu +1 - pmu;
%     plot(F_rescaled(:,mu), '--','Color', scale_gradient(pmu+1,:))
    plot(Fnew(:,mu), 'Color', scale_gradient(mu,:))
end

%% Get maxes of new
for mu = 1:nMu
[Fnewmax(:,mu), indnewmax(:,mu)] = get_max(Fnew(:,mu), window_len/dt);
end

%% Plotting
nppts = 100;
figure(7)
clf
subplot(2,2,1)
hold on
pwindow = time < window_len * 2;
sset = ceil(linspace(find(pwindow,1,'first'),find(pwindow,1,'last'),nppts));
plot(time(sset)/1000,sum(Fnew(sset,:),2), 'o', 'color', bloo, 'markersize',5,'linewidth',1.1)
plot(time(pwindow)/1000,F_Ding03(pwindow), '-', 'color', reed, 'linewidth', 1.5)
% plot(time(indnewmax(1:2,1))/1000, sum(Fnewmax(1:2,:),2),'^', 'color', bloo)
% plot(time(imaxD03(1:2))/1000, fmaxD03(1:2),'o', 'color', reed)
text(0.05,1,char('A'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 20)
title('First Two Trains')
ylabel('Force (N)')
xlabel('Time (s)')
legend(["Motor Unit Model" "Ding Model"])
axis([0 inf 0 max(F_Ding03)*1.5])

subplot(2,2,2)
hold on
pwindow = time > window_len * (pulse_trains - 2);
sset = ceil(linspace(find(pwindow,1,'first'),find(pwindow,1,'last'),nppts));
plot(time(sset)/1000,sum(Fnew(sset,:),2), 'o', 'color', bloo, 'markersize', 5,'linewidth',1.1)
plot(time(pwindow)/1000,F_Ding03(pwindow), '-', 'color', reed, 'linewidth', 1.5)

% plot(time(indnewmax(pulse_trains-1:end,1))/1000, sum(Fnewmax(pulse_trains-1:end,:),2),'^', 'color', bloo)
% plot(time(imaxD03(pulse_trains-1:end))/1000, fmaxD03(pulse_trains-1:end),'o', 'color', reed)
text(0.05,1,char('B'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 20)
title('Last Two Trains')
ylabel('Force (N)')
xlabel('Time (s)')
legend(["Motor Unit Model" "Ding Model"])
axis([-inf inf 0 max(F_Ding03)*1.5])

subplot(2,2,3)
hold on
plot(time(indnewmax(:,1))/1000, sum(Fnewmax,2),'o', 'color', bloo,'linewidth',1.1, 'markersize', 4)
plot(time(imaxD03)/1000, fmaxD03,'-', 'color', reed,'linewidth',1.5)
text(0.05,1,char('C'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 20)
title('Maximum Force -- Original vs Parameterized')
ylabel('Force (N)')
xlabel('Time (s)')
axis([-inf inf 0 max(F_Ding03)*1.5])
legend(["Motor Unit Model" "Ding Model"])

subplot(2,2,4)
hold on
plot(time(imaxD03)/1000, fmaxD03 - sum(Fnewmax,2),'.', 'color', reed,'markersize', 9)
line([0 max(time(imaxD03)/1000)], [0 0 ], 'linestyle', '--')
text(0.05,1,char('D'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 20)
text(0.17,0.79,char('Motor Unit Model < Ding Model'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 11)
text(0.17,0.28,char('Motor Unit Model > Ding Model'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 11)

title('Maximum Force Residuals')
ylabel('Maximum Force Residual (N)')
xlabel('Time (s)')
axis([-inf inf -40 40])
% subplot(2,2,4)
% hold on
% for mu = 1:nMu
%     pmu = nMu - mu +1
%     plot(time(indnewmax(:,1))/1000, force_prop(:,pmu), 'color', scale_gradient(pmu,:))
%     hold on
% end
% text(0.05,1,char('B'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 20)
% title('Maximum Force -- Original vs Parameterized')
% ylabel('Force (N)')
% xlabel('Time (s)')
% axis([-inf inf 0 max(F_Ding03)*1.5])

sizemod = 1.6;
set(gcf, 'Position', [0, 0, 600*sizemod, 400*sizemod])
export_fig 'C:\Users\u153094\Documents\UPF\GradSchool\Thesis Project\Paper Files\Figures\ParameterizedvsOriginalForce.png' -transparent -m3

%% Plot % decline and summed forces
figure(8)
clf
subplot(2,1,1)
pwindow = time < window_len * 2;
for mu = 1:nMu
    plot(time(pwindow)/1000, Fnew(pwindow,mu), '-','color', scale_gradient(mu,:),'linewidth',1.1)
    hold on
end
plot(time(pwindow)/1000,sum(Fnew(pwindow,:),2), '--', 'color', bloo,'linewidth',1.1)
% plot(time(pwindow)/1000,F_Ding03(pwindow), '--', 'color', reed)
legend(['MU' + string([1:nMu]), "Muscle"],'Location', 'eastoutside')
axis([-inf inf 0 max(F_Ding03)*1.1])
text(0.015,0.95,char('A'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 20)
title('Motor Unit and Muscle Forces')
ylabel('Force (N)')
xlabel('Time (s)')

subplot(2,1,2)
for mu = 1:nMu
    plot(time(indnewmax(:,1))/1000, force_prop(:,mu)*100, '-','color', scale_gradient(mu,:),'linewidth',1.1)
    hold on
end
plot(time(indnewmax(:,1))/1000, percent_decline*100, '--','color', bloo,'linewidth',1.1)
legend(['MU' + string([1:nMu]), "Muscle"],'Location', 'eastoutside')
axis([-inf inf 0 100])
text(0.015,.25,char('B'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 20)
title('Motor Unit and Muscle Fatigue Curves')
ylabel('Percent of Initial Maximum Force')
xlabel('Time (s)')

sizemod = 1.6;
set(gcf, 'Position', [0, 0, 600*sizemod, 400*sizemod])
% export_fig 'C:\Users\u153094\Documents\UPF\GradSchool\Thesis Project\Paper Files\Figures\PercentDecline.png' -transparent -m3
