function [time, stim_times, pulses, train_starts, train_dur] = create_pulsetrain(pulse_trains, nMu, dt, ptype, stype)

if nargin == 1
    nMu = 1;
    dt = 0.1;
    ptype = 'v50c100';
    stype = [];
elseif nargin == 2
    dt = 0.1;
    ptype = 'v50c100';
    stype = [];
elseif nargin == 3
    ptype = 'v50c100';
    stype = [];
elseif nargin == 4   
    stype = [];
end

pulsewidth = 1; %ms
doublet_spacing = 5; %ms

if strcmp(ptype, 'v50c100')
    v50pulses = 6;
    v_ipi = 50;
    c100pulses = 6;
    c_ipi = 100;
    delay = 500; %ms
    
    % set n first pulses
    v50stim_times = (0: v_ipi : v_ipi * (v50pulses - 1)) ;
    
    % add doublet
    v50stim_times = [v50stim_times(1), v50stim_times(1) + pulsewidth + doublet_spacing, v50stim_times(2:end)];
    
    c100pulse_idxs = (0: c_ipi : c_ipi * (c100pulses-1)) + delay + v50stim_times(end);
    stim_times = [v50stim_times , c100pulse_idxs];
    train_dur = stim_times(end) + delay;
    
    tstim_times = stim_times;
    for irep = 1:pulse_trains-1
        tstim_times = [tstim_times , stim_times + train_dur * irep];
    end
    
    stim_times = tstim_times';
    
    tfinal = train_dur * (pulse_trains); %Final time (ms)
    %     time(pulse_idxs(end)) - time(pulse_idxs(1))
    train_starts = 0:train_dur:train_dur * (pulse_trains-1);
    
    pulses = length(stim_times);
    time = 0:dt:tfinal; %time scale (ms)
    
    %Make identical pulse trains for each motor unit
    for mu = 2:nMu
        stim_times(:,mu) = stim_times(:,1);
    end
elseif strcmp(ptype, 'c100')
    train_dur = 100;
    stim_times = 0:100:100*(pulse_trains-1);
    tfinal = train_dur * (pulse_trains+5);
    train_starts = 0:train_dur:train_dur * (pulse_trains-1);
    pulses = length(stim_times);
    time = 0:dt:tfinal;
    for mu = 2:nMu
        stim_times(:,mu) = stim_times(:,1);
    end
elseif strcmp(ptype, 'c33')
    train_dur = 2000;
    ipi = 33;
    stim_times = 0:ipi:1000;
    tstim_times = stim_times;
    for irep = 1:pulse_trains-1
        tstim_times = [tstim_times , stim_times + train_dur * irep];
    end
    stim_times = tstim_times';
    tfinal = train_dur * (pulse_trains);
    train_starts = 0:train_dur:train_dur * (pulse_trains-1);
    pulses = length(stim_times);
    time = 0:dt:tfinal;
    for mu = 2:nMu
        stim_times(:,mu) = stim_times(:,1) + (mu -1)*ipi ;
    end
end
time = time(1:end-1);

if strcmp(stype, 'interleaved')
    for mu = 1:nMu
    stim_times(:,mu) = stim_times(:,mu) + 2000*mu/(nMu+1);%floor(ipi*mu/(nMu+1));
    end
end