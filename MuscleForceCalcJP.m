%% Pulsetrain Definition
dt = 0.1; %ms
nMu = 10; %number motor units
pulse_trains = 75;
window_len = 1750;
mode = 'v50c100';

[time, stim_times, pulses, train_starts] = create_pulsetrain(pulse_trains, nMu, dt);
LoadColorScales;
global scale_gradient
global scale_discrete
% train_dur = train_starts(2);
%% Force Calculations -- Setup

%Time to peak (P) and max twitch force(T) from Fugelvand 1993  
forceRange = 120; % 10-fold force diffecence between largest and smallest MU
fatigueRange = 180; % range of fatigue rates across the motor units (300 best)
contractionRange = 3; % Range in contraction times,
ct_max = 90;
% FatFac = 0.0225;    % fatigue factor (FF/S) percent of peak force of MU per second 
tick_breaks = 3;
[P, fatigue, arborLengths, m_csa, ct]  = define_motor_pool(nMu, forceRange, fatigueRange, contractionRange);
fatigue = fatigue / min(fatigue);

figure(1)
clf
% Plot Peak twitch force
subplot(3,2,1)
hold on
box off
plot(P, '.', 'MarkerSize', 8, 'Color', bloo)
text(0.05,0.95,char('A'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 20)
xlabel('Motor Unit Number')
ylabel('Normalized Twitch Force')
title('Twitch Force vs Motor Unit')
xticks([0:nMu/tick_breaks:nMu])
yticks([0:nMu/tick_breaks:nMu])
axis([0 nMu 0 forceRange])

% Hist Force Values
subplot(3,2,2)
histogram(P, 10, 'Normalization','probability', 'FaceColor', bloo)
set(gca,'GridLineStyle','none')
text(0.05,0.95,char('C'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 20)
xlabel('Normalized Twitch Force')
ylabel('Proportion of Motor Units')
title('Motor Unit Force Distribution')
box off
xticks([0:nMu/tick_breaks:nMu])
axis([0 inf 0 1])

% Fatigue Factor vs Twitch Force
subplot(3,2,3)
plot(P, fatigue, '.', 'MarkerSize', 8, 'Color', bloo)
text(0.05,0.95,char('B'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 20)
xlabel('Normalized Twitch Force')
ylabel('Normalized Fatigue Factor')
title('Fatigue Factor vs Twitch Force')
box off
xticks([0:nMu/4:nMu])
yticks([0:fatigueRange/tick_breaks:fatigueRange])
axis([0 120 0 fatigueRange])

%Contraction time vs Twitch force
subplot(3,2,4)
hold on
box off
plot(P', ct, '.', 'MarkerSize', 8, 'Color', bloo)
text(0.05,0.95,char('D'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 20)
xlabel('Normalized Twitch Force')
ylabel('Contraction Time (ms)')
title('Contraction Time vs Twitch Force')
xticks([0:nMu/tick_breaks:nMu])
yticks([0:ct_max/tick_breaks:ct_max])

axis([0 120 0 inf])

% Arbor Length vs Twitch force 
subplot(3,2,5)
hold on
box off
plot(P', arborLengths, '.', 'MarkerSize', 8, 'Color', bloo)
text(0.05,0.95,char('E'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 20)
xlabel('Normalized Twitch Force')
ylabel('Normalized Arbor Length')
title('Arbor Length vs Twitch Force')
xticks([0:nMu/tick_breaks:nMu])
axis([0 120 0 1])


% CSA vs Twitch force 
subplot(3,2,6)
hold on
box off
plot(P, m_csa, '.', 'MarkerSize', 8)
text(0.05,0.95,char('F'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 20)
xlabel('Normalized Twitch Force')
ylabel('Normalized Axonal Diameter')
title('Axonal Diameter vs Twitch Force')
xticks([0:nMu/tick_breaks:nMu])
axis([0 120 0 1])

sizemod = 1.15;
set(gcf, 'Position', [400, 0, 600*sizemod, 600*sizemod])
% export_fig 'C:\Users\u153094\Documents\UPF\GradSchool\Thesis Project\Paper Files\Figures\MotorPool6.png' -transparent -m3
%% Force Calculations --  Plot Optimized Doll Values
param_names = ["A_rest" "t_fat" "t1" "t2" "aA" "aT1" "aKm"];
param_vals = [1.44; 139.74*1000; 25.92; 84.27; -4.15e-7; 2e-5; 3e-8];
[F, time, Cn, Km, R0, A, T1] = calc_MU_force(1, 75, dt, param_names, param_vals);

%% Force Calculations --  Plot Original Ding 2002 Values
param_names = ["A_rest" "t_fat" "t1" "t2" "aA" "aT1" "aKm"];
param_vals = [5.1; 53.4*1000; 43.8; 124.4; -8.8e-7; 5.7e-6; 1.9e-8];
pulse_trains = 75;
window_len = 1750;
[F, time, Cn, Km, R0, A, T1] = calc_MU_force(1, 75, dt, param_names, param_vals);

%% Force Calculations --  Plot 2003 Ding Values
param_names = ["A_rest" "t_fat" "t1" "t2" "aA" "aT1" "aKm"];
t12scale = .75;
pulse_trains = 120;
param_vals = [3.009; 127*1000; 43.8*t12scale; 124.4*t12scale; -4.15e-7; 5.7e-6; 1.9e-8];
[F, time, Cn, Km, R0, A, T1] = calc_MU_force(1, pulse_trains, dt, param_names, param_vals, 'c33');
clf
window_len = 2000;
plot(time/1000,F)
%% Plot reconstruction
figure(2)
clf
hold on

% subplot(2,6,[1 2])
% plot(time/1000,A)
% ylabel('A (N/ms)')
% xlabel('Time (s)')
% axis([0 inf 0 max(A)*1.2])
% 
% subplot(2,6,[3 4])
% plot(time/1000,R0)
% ylabel('Km (Unitless)')
% xlabel('Time (s)')
% axis([0 inf 0 max(R0)*1.2])
% 
% subplot(2,6,[5 6])
% plot(time/1000,T1)
% ylabel('T1 (ms)')
% xlabel('Time (s)')
% axis([0 inf 0 max(T1)*1.2])
[pmax, ipmax] = get_max(F, window_len/dt);

subplot(1,3,1)
hold on
pwindow = time < window_len * 2;
plot(time(pwindow)/1000,F(pwindow))
plot(time(ipmax(1:2))/1000, pmax(1:2),'o', 'color', reed)
text(0.05,1,char('A'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 20)
title('First Two Trains')
ylabel('Force (N)')
xlabel('Time (s)')
axis([0 inf 0 max(F)*1.2])

subplot(1,3,2)
hold on
pwindow = time > window_len * (pulse_trains - 2);
plot(time(pwindow)/1000,F(pwindow))
plot(time(ipmax(pulse_trains-1:end))/1000, pmax(pulse_trains-1:end),'o', 'color', reed)
text(0.05,1,char('B'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 20)
title('Last Pulse Trains')
ylabel('Force (N)')
xlabel('Time (s)')
axis([-inf inf 0 max(F)*1.2])

subplot(1,3,3)
hold on
plot(time(ipmax)/1000, pmax,'.', 'color', reed)
text(0.05,1,char('C'),'Units', 'Normalized', 'VerticalAlignment', 'Top', 'FontSize', 20)
title('Peak Force vs Time')
ylabel('Peak Force (N)')
xlabel('Time (s)')
xticks([0:25:max(time)/1000])
axis([0 inf 0 max(F)*1.2])

sizemod = 1.15;
set(gcf, 'Position', [0, 0, 1000*sizemod, 200*sizemod])
% export_fig 'C:\Users\u153094\Documents\UPF\GradSchool\Thesis Project\Paper Files\Figures\Ding2002Force.png' -transparent -m3

%% Plotting
% csvwrite('Forces', horzcat(time', Fsum(2:end), F(2:end,:)))
% scatter3(gcolors(:,1), gcolors(:,2), gcolors(:,3.),30, gcolors/255)
window_length = 1750;
figure(2)
clear h
clf
hold on

for pmu = 1:nMu
    mu = nMu + 1 - pmu;
    %Show instantaneous force curves, comment out for max force curves
    h(mu*2 -1) = plot(time/1000,F(:,mu), 'Color', scale_gradient(mu,:));
    %Calculate Maximum force in each train
    [fmax(:,mu), indmax(:,mu)] = get_max(F(:,mu), window_length/dt);
    %plot a dot at each max
    h(mu*2) = plot(indmax(:,mu)/1000*dt, fmax(:,mu), 'o', 'Color',scale_gradient(mu,:));
end
% gridxy(time(stim_times(:,1))
title('Individual Motor Unit Force vs Time')
xlabel('Time (s)')
ylabel('Force (N)')
legend(h([1:2:nMu*2]), ['MU' + string([1:nMu])])
% legend(h([1,11:2*10:nMu*2]), ['MU' + string([1, 10:10:nMu])])
hold off
% axis([0 4 0 15])
clear h

% figure(3)
% clf
% hold on
% plot(time/1000,Fsum(2:end),'Color',[0,0,0])
% for mu = 1 :nMu
%     %Show instantaneous force curves, comment out for max force curves
%     h(mu) = plot(time/1000,F(2:end,mu), 'Color', scale_gradient(mu,:));
% end
% % gridxy(time(pulse_idxs))
% legend(['Summed Force','MU' + string(1:nMu)])
% title('Summed and Motor Unit Forces')
% xlabel('Time (s)')
% ylabel('Force (N)')
% hold on

% figure(3)
% clf
% hold on
% subplot(1,3,1)
% plot(time/1000,A(2:end,1),'.')
% axis([0,inf,0,ceil(max(A(2:end,1))*1.1)])
% % gridxy(time(pulse_idxs))
% title('A Parameter Fatigue')
% xlabel('Time (s)')
% ylabel('A')
% 
% subplot(1,3,2)
% plot(time/1000,R0(2:end,1),'.')
% axis([0,inf,0,ceil(max(R0(2:end,1))*1.1)])
% title('R0 Parameter Fatigue')
% xlabel('Time (s)')
% ylabel('R0')
% 
% subplot(1,3,3)
% plot(time/1000,T1(2:end,1),'.')
% axis([0,inf,0,ceil(max(T1(2:end,1))*1.1)])
% title('T1 Parameter Fatigue')
% xlabel('Time (s)')
% ylabel('T1')

percent_decrease = 1 - fmax(end,:) ./ fmax(1,:)
percent_decrease_per_second = (1 - fmax(end,:) ./ fmax(1,:)) / (time(end)/1000);

%% Compare multi to single
Fsingle = csvread('F1MU75Trains.csv');
figure(4)
clf
hold on 

plot(time,Fsingle(:,1), 'Color', bloo)
plot(time,Fsum(2:end), '--','Color', reed)
d_mse(Fsingle(:,1), Fsum)
legend('Single Motor Unit', string(nMu) + ' Motor Units')
% Digitized Values : 1st max = 326.143791, last max = 169.281046