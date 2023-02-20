function [output] = granule_STDP(net_params,sim_params,options)
%granule_STDP Generates a population of synapses with a kesten process, multiplicative STDP, and heterosynaptic scaling. 

% Simulation parameters define simulation time
if nargin<2 || isempty(sim_params)
  sim_params=struct;
end
if ~isfield(sim_params,'dt') % Timestep size (10ms)
    sim_params.dt=0.0005;
end
if ~isfield(sim_params,'max_time') % Maximum time to receive inputs
    sim_params.max_time=200;
end
if ~isfield(sim_params,'silent') % Just silent part?
    sim_params.silent=false;
end
if ~isfield(sim_params,'record_filter') % How often to record weights?
    sim_params.record_filter={false,10000};
end

% Network parameters define network structure
if nargin<1 || isempty(net_params)
  net_params=struct;
end
if ~isfield(net_params,'n_neuron') % Number of neurons (repetitions)
    net_params.n_neuron=500;
end
if ~isfield(net_params,'n_spine')  % Number of spines per neuron
    net_params.n_spine=100;
end
if ~isfield(net_params,'tot_size') % Total spine volume/EPSP size (Effectively in mV)
    net_params.tot_size=100; 
end
if ~isfield(net_params,'pot_rate') % Rate of potentiation with STDP
    net_params.pot_rate=0.000035; 
end
if ~isfield(net_params,'dep_rate') % Relative rate of depression with STDP
    net_params.dep_rate=0.000000015; 
end
if ~isfield(net_params,'tau_STDP') % Timescale of STDP
    net_params.tau_STDP=0.5; 
end
if ~isfield(net_params,'tau_w') % Timescale of adaptation
    net_params.tau_w=0.1; 
end
if ~isfield(net_params,'w') % Strength of adaptation
    net_params.w=10; 
end

if ~isfield(net_params,'heterosynaptic') % Does heterosynaptic plasticity keep total weight constant?
    net_params.heterosynaptic=true; 
end
if ~isfield(net_params,'weight_adaptation') % Does the threshold adjust to maintain spike rates? What rate?
    net_params.weight_adaptation={false,100}; 
end
if ~isfield(net_params,'kesten_noise') % Strength of multiplicative (Kesten) noise process (in %)
    net_params.kesten_noise=2*sqrt(sim_params.dt)/0.3162; 
end
if ~isfield(net_params,'gauss_noise') % Strength of additive (gaussian) noise process (in mV)
    net_params.gauss_noise=0*sqrt(sim_params.dt)/0.3162; 
end
if ~isfield(net_params,'threshold') % Firing threshold (mV)
    net_params.threshold=10;
end
if ~isfield(net_params,'time_constant') % Membrane time constant (s)
    net_params.time_constant=0.01;
end
if ~isfield(net_params,'input_pattern') % Are inputs Poisson or periodic?
    net_params.input_pattern='Poisson'; % 'periodic' for regular spiking
end
if ~isfield(net_params,'input_rate') % Rate of inputs per feedforward cell per timestep
    net_params.input_rate=200;
end
if ~isfield(net_params,'input_corr') % Cross-correlation between inputs (between 0 for independent and 1 for fully synchronous)
    net_params.input_corr=0;
end
if ~isfield(net_params,'partial_stim') % Fraction of synapses that receive inputs, set below 1 for Tassilo's data...
    net_params.partial_stim=1;
end
if ~isfield(net_params,'init_dist') % Is there an initial distribution of weights?
    net_params.init_dist=[]; % No, other options are a vector of size (n_neuron*n_spine,1) or (n_spine,n_neuron);
end
if ~isfield(sim_params,'pre_run')
    if isempty(net_params.init_dist) % Time before input
        sim_params.pre_run=50;
    else
        sim_params.pre_run=0;
    end
end


% Options (plot distributions)
if nargin<3 || isempty(options)
  options=''; % '-S1' plots initial and final distributions, '-S2' plots all intermediate
end

% Pre-run spine sizes
if  isempty(net_params.init_dist)  % Uniform start (default)
    spine_sizes=net_params.tot_size/net_params.n_spine*ones(net_params.n_spine,net_params.n_neuron);
elseif size(net_params.init_dist,2)==1  % Single vector for different distributions
    raw_sizes=net_params.init_dist(randperm(length(net_params.init_dist))); % Reshuffle inputs
    spine_sizes=reshape(raw_sizes,net_params.n_spine,net_params.n_neuron);
    spine_sizes=net_params.tot_size*spines_sizes./sum(spine_sizes);
elseif size(net_params.init_dist,2)==net_params.n_neuron % Using final distribution of other parameters
    spine_sizes=net_params.init_dist;
end

% Silent period
if sim_params.pre_run>0
    for dt_ind=1:round(sim_params.pre_run/sim_params.dt)
        spine_sizes=spine_sizes.*normrnd(1,net_params.kesten_noise/100,net_params.n_spine,net_params.n_neuron)+normrnd(0,net_params.gauss_noise,net_params.n_spine,net_params.n_neuron); % Stochastic growth
        if net_params.heterosynaptic
            spine_sizes=net_params.tot_size*spine_sizes./sum(spine_sizes); % Heterosynaptic part
        end
    end
end

grid_len=ceil((sim_params.max_time+sim_params.dt/2)/sim_params.dt);
output=struct;
output.initial_dist=spine_sizes;
nSP=zeros(grid_len,1); % Record spikes at each timestep
if sim_params.record_filter{1} % Record intermediate weights
    weights_output=zeros(net_params.n_spine*net_params.n_neuron,round(sim_params.max_time/(sim_params.dt*sim_params.record_filter{2})));
end

if ~sim_params.silent
    % Establish connectivity
    n_each=round(net_params.partial_stim*net_params.n_spine); % How many synapses receive input?

    % Generate and regularise spikes
    if strcmp(net_params.input_pattern,'Poisson') % Poisson inputs
        input_spiketimes=Poisson_correlated(sim_params.max_time, n_each,net_params.input_rate,net_params.input_corr);
    elseif strcmp(net_params.input_pattern,'periodic') % Periodic inputs
        input_spiketimes=periodic_correlated(sim_params.max_time, n_each,net_params.input_rate,net_params.input_corr);
    else
        error('Invalid stimulation pattern');
    end
   
    binned_spikes=zeros(net_params.n_spine,grid_len);
    for ward=1:n_each
        nspike=length(input_spiketimes{ward});
        for sooth=1:nspike
            binned_spikes(ward,ceil(input_spiketimes{ward}(sooth)/sim_params.dt))=1;
        end
    end

    net_E=net_params.threshold*rand(net_params.n_neuron,1); % Neuronal activities (initially uniform between reset and threshold)
    pre_times=zeros(net_params.n_spine,1); % Time since last EPSP at each spine
    post_times=zeros(net_params.n_neuron,1); % Time since last spike at each neuron
    for dt_ind=1:round(sim_params.max_time/sim_params.dt)
        net_E=net_E*exp(-sim_params.dt/net_params.time_constant); % Mebrane leak
        % Inputs
        net_E=net_E+spine_sizes'*binned_spikes(:,dt_ind)-net_params.w*exp(-pre_times'/net_params.tau_w);
        pre_times(binned_spikes(:,dt_ind)==1)=0;
        pre_times=pre_times+sim_params.dt*ones(net_params.n_spine,1); % Increment times since EPSPs

        % Record spikes
        post_activity=net_E>=net_params.threshold; % E cells above threshold
        net_E(post_activity)=0; % Reset after spikes
        %
        nSP(dt_ind)=nnz(post_activity); % Record rate

        % STDP
        if nnz(post_activity)>0 % Check for spiking
            pre_tiled=repmat(exp(-pre_times/net_params.tau_STDP),1,net_params.n_neuron);
            spine_sizes=spine_sizes+net_params.pot_rate*pre_tiled; % Potentiate
        end
        if nnz(binned_spikes(:,dt_ind))>0 % Check for EPSPs
            post_tiled=repmat(exp(-post_times'/net_params.tau_STDP),net_params.n_spine,1);
            spine_sizes=spine_sizes.*(ones(net_params.n_spine,net_params.n_neuron)-net_params.dep_rate*post_tiled); % Depress
        end
        post_times=post_times+sim_params.dt*ones(net_params.n_neuron,1); % Increment times since spikes

        % Intrinsic noise
        spine_sizes=spine_sizes.*normrnd(1,net_params.kesten_noise/100,net_params.n_spine,net_params.n_neuron)+normrnd(0,net_params.gauss_noise,net_params.n_spine,net_params.n_neuron); % Stochastic growth
        if net_params.heterosynaptic
            spine_sizes=net_params.tot_size*spine_sizes./sum(spine_sizes); % Heterosynaptic part
        end
        if sim_params.record_filter{1} % Record intermediate weights
            if rem(dt_ind,sim_params.record_filter{2})==0
                weights_output(:,round(dt_ind/sim_params.record_filter{2}))=spine_sizes(:);
            end
        end
    end
    output.fin_dist=spine_sizes;
    output.spike_rates=nSP;

    if contains(options,'-S1') % Plot inital and final distributions
        [ni,Ci]=hist(output.initial_dist(:),100);
        ni=ni/trapz(Ci,ni); % Normalise
        [nf,Cf]=hist(output.fin_dist(:),100);
        nf=nf/trapz(Cf,nf); % Normalise
        figure
        hold on
        plot(Ci,ni,'black')
        plot(Cf,nf,'blue')
    end
    if contains(options,'-S2') % Plot intermediate distributions if recorded (see sim_params.record_filter)
        cscale=[linspace(0.7,0,size(weights_output,2)) ; linspace(0.7,0,size(weights_output,2)) ; linspace(0.7,0,size(weights_output,2))];
        figure
        hold on
        for run_ind=1:size(weights_output,2)
            [n,C]=hist(weights_output(:,run_ind),100);
            n=n/trapz(C,n); % Normalise
            plot(C,n,'Color',cscale(:,run_ind)')
        end
    end
end

end


function [input_spiketimes] = Poisson_correlated(max_time,n_spine,rate,corr)
%POISSON_CORRELATED Generates a set of Poisson correlated spike times with correlations between trains.

input_spiketimes=cell(n_spine,1);
if corr>0
    master_rate=rate/corr;    
    master_spiketimes=zeros(1);
    running_time=exprnd(1/master_rate);
    n_spike=0;
    while running_time<max_time
        n_spike=n_spike+1;
        master_spiketimes(n_spike)=running_time;
        running_time=running_time+exprnd(1/master_rate);
    end
    for spine_ind=1:n_spine
        these_spiketimes=zeros(1);
        spikes_here=0;
        for spike_ind=1:n_spike
            tst=rand(1);
            if tst<corr
                spikes_here=spikes_here+1;
                these_spiketimes(spikes_here)=master_spiketimes(spike_ind)+0.005*randn(1);
            end
        end
        input_spiketimes{spine_ind}=these_spiketimes;
    end
else
    for spine_ind=1:n_spine
        these_spiketimes=zeros(1);
        spikes_here=0;
        running_time=exprnd(1/rate);
        while running_time<max_time
            spikes_here=spikes_here+1;
            these_spiketimes(spikes_here)=running_time;
            running_time=running_time+exprnd(1/rate);
        end
        input_spiketimes{spine_ind}=these_spiketimes;
    end
end
end

function [input_spiketimes] = periodic_correlated(max_time,n_spine,rate,corr)
%PERIODIC_CORRELATED Generates a set of periodic correlated spike times with correlations between trains.

input_spiketimes=cell(n_spine,1);
if corr>0
    master_rate=rate/corr;
    master_spiketimes=zeros(1);
    running_time=1/master_rate;
    n_spike=0;
    while running_time<max_time
        n_spike=n_spike+1;
        master_spiketimes(n_spike)=running_time;
        running_time=running_time+1/master_rate;
    end

    for spine_ind=1:n_spine
        these_spiketimes=zeros(1);
        spikes_here=0;
        for spike_ind=1:n_spike
            tst=rand(1);
            if tst<corr
                spikes_here=spikes_here+1;
                these_spiketimes(spikes_here)=master_spiketimes(spike_ind)+0.00005*randn(1);
            end
        end
        input_spiketimes{spine_ind}=these_spiketimes;
    end
else
    for spine_ind=1:n_spine
        these_spiketimes=zeros(1);
        spikes_here=0;
        running_time=exprnd(1/rate); % Shuffle starts to break correlation
        while running_time<max_time
            spikes_here=spikes_here+1;
            these_spiketimes(spikes_here)=running_time;
            running_time=running_time+1/rate;
        end
        input_spiketimes{spine_ind}=these_spiketimes;
    end
end
end