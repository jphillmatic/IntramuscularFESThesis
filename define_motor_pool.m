function [P, t_fat, arborLengths, m_csa, ct]  = define_motor_pool(nMu, forceRange, fatigueRange, contractionRange)
   
n = 1:1:nMu;

if nMu == 1
    t_fat = 1;
    P = 1;
    m_csa = 1;
    arborLengths = 1;
    ct = 68.8;
    return
end

%Time to peak (P) and max twitch force(T) from Fugelvand 1993
if nargin == 1
    forceRange = 120; % 10-fold force diffecence between largest and smallest MU
    fatigueRange = 180;% range of fatigue rates across the motor units (300 best)
    contractionRange = 3; % range of time to peack twitch force
end

%Calculate Force Values
P = exp( log(forceRange).* (n-1) / (nMu-1)); %Unitless

tL = 145; % longest contraction time (90ms)

c = log(forceRange)/log(contractionRange);% scale factor
ct = tL * (1./P).^(1/c); %contraction time of motor units (ms)


% Fatigue Rate from paper
FatFac = 0.000125;
t_fat = FatFac * exp(log(fatigueRange)*(n-1)/(nMu-1));

%Cross sectional area
% TODO FIND CSA OF A MOTOR UNIT AND SCALE
m_csa = P.^(0.4544);
m_csa = m_csa / sum(m_csa)*3.05*10e-3;

% Arbor Length
% TODO -- FIND ARBOR LENGTH OF A MOTOR UNIT AND SCALE
arborLengths = P.^(0.4938);
arborLengths = arborLengths / max(arborLengths);

% Potvin fatigue 2017+, updated but unpublished?
% Fuglevand defined fatigue rate for largest motor unit
% FatFac = 0.0225;    % fatigue factor (FF/S) percent of peak force of MU per second 
% fatigue rate for each motor unit FROM FUGELVAND CODE
% b2 = log(fatigueRange)/(nMu-1);
% mufatrate = exp(b2 .* (n-1));
% fatigue = mufatrate * (FatFac / fatigueRange) .* P;