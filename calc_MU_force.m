function [Force, time, Cn, Km, R0, A, T1, train_dur] = calc_MU_force(nMu, pulse_trains, dt, param_names, muscle_params, stim_type, p_type)
%% Force Calculations -- Setup

if nargin <= 3
    muscle_params  = [];
    param_names = [];
    stim_type = "v50c100";
    p_type = [];
elseif nargin <= 5
    stim_type = "v50c100";
    p_type = [];
elseif nargin <= 6
    p_type = [];
end

[time, stim_times, pulses, train_dur] = create_pulsetrain(pulse_trains, nMu, dt, stim_type, []);

% stim_times = timbos;

%Time to peak (P) and max twitch force(T) from Fugelvand 1993
forceRange = 120; % 10-fold force diffecence between largest and smallest MU
fatigueRange = 180; % range of fatigue rates across the motor units (300 best)
contractionRange = 3; % Range in contraction times,
ct_center = 68.8; %ms to peak contraction time of original muscle model

%Define twitch force distribution
P = define_motor_pool(nMu, forceRange, fatigueRange, contractionRange);

%Stationary Parameters
t1 = 43.8; %Time constant of force decline in absence of actin linkages (ms)
%<- nonfatigued 43.8, fatigued = 89.2
t2 = 124.4; %Time constant of force decline with strong actin linkages (ms)
%<- nonfatigued, fatigued = 1564.5
%paralyzed; 58-78ms, nonparalyzed; 124-1564 ms


A_total = 5.1;
A_rest = P / sum(P) * A_total; %A at rest

t_fat_center = 53.4*1000;
t_fat = repmat(t_fat_center,1,nMu); %(s -> ms)Fatigue Rate 47.8, 46.9, 53.4, 126 from lit

% R0 = repmat(1.3,1,nMu); %R0 at rest
tc = 20; %tc at rest
Km_rest =  ones(1,nMu)*0.3;

aA = -8.8e-7 * ones(1,nMu);
aT1 = 5.7e-6 * ones(1,nMu);
aKm = 1.9e-8 * ones(1,nMu);

if ~isempty(muscle_params)
    for pnm = 1:length(param_names)
        switch param_names(pnm)
            case "A_rest"
                A_rest = muscle_params(pnm,:);
            case "t_fat"
                t_fat = muscle_params(pnm,:);
            case "t1"
                t1 = muscle_params(pnm,:);
            case "t2"
                t2 = muscle_params(pnm,:);
            case "aA"
                aA = muscle_params(pnm,:);
            case "aT1"
                aT1 = muscle_params(pnm,:);
            case "aKm"
                aKm = muscle_params(pnm,:);
            case "tc"
                tc = muscle_params(pnm,:);
            case "Km_rest"
                Km_rest =  muscle_params(pnm,:);
        end
    end
end

%% Force Calculations -- Loop
% Fsum = zeros(length(time),1);

%Tracked Values
Cn = zeros(length(time),nMu); %Concentration of calcium
Force = zeros(length(time),nMu); % force generated

A = zeros(length(time), nMu);
A(1,:) = A_rest; %Scaling Factor for muscle force

Km = zeros(length(time), nMu);  %Gradient flow constant
Km(1,:) = Km_rest;

T1 = zeros(length(time), nMu);  %Gradient flow constant
T1(1,:) = ones(1,nMu).*t1;

R0 = zeros(length(time), nMu); % Force scaling from potentiation

%define time between each pulse for R calculations
pulse_space = stim_times - wshift('2D',stim_times,[-1,0]);
pulse_space(1,:) = 10000;

for i = 1:length(time)-1
    
    R0(i,:) = Km(i,:) + 1.04;
    
    Rimp = zeros(pulses,nMu); %Strength increase due to repeated stimulation

    
%     for mu = 1:nMu
%         initiatedPulses = (time(i) >= stim_times);
%         lp = max(sum(initiatedPulses,1));
%         Rimp(mu) = (1 + (R0(i,:) - 1) .* exp( -pulse_space/ tc )) .* exp( - (time(i) - stim_times) /tc);
%     end
    
    initiatedPulses = (time(i) >= stim_times);
    lp = max(sum(initiatedPulses,1));
    Rimp = (1 + (R0(i,:) - 1) .* exp( -pulse_space/ tc )) .* ... % Ri
        exp( - (time(i) - stim_times) /tc);
    
    %Calcium calculation
    dCn = sum(Rimp(1:lp,:),1) / tc - Cn(i,:) / tc;
    
    %Force calculation
    dF = A(i,:) .* Cn(i,:) ./ (Km(i,:) + Cn(i,:)) - Force(i,:) ./ (T1(i,:) + t2 .* (Cn(i,:) ./ (Km(i,:) + Cn(i,:))));
    
    %Calculate change in fatiguing parameters
    dT1 = (- ((T1(i,:) - t1) ./ t_fat) + aT1 .* Force(i,:))*dt;
    dA = (- ((A(i,:) - A_rest) ./ t_fat) + aA .* Force(i,:))*dt;
    dKm = (- ((Km(i,:) - Km_rest) ./ t_fat) + aKm .* Force(i,:)) *dt;
    
    %Set values for next iter
    Km(i+1,:) = Km(i,:) + dKm;
    A(i+1,:) = A(i,:) + dA;
    T1(i+1,:) = T1(i,:) + dT1;
    
    Cn(i+1,:) = Cn(i,:) + dCn * dt;
    Force(i+1,:) = Force(i,:) + dF * dt;
    
        
end
% Fsum = sum(Force,2);