% Spiking neuron model for GLM identification

clear all; close all
addpath(genpath('/Users/yanj11/code/Shotgun'))

% dimension
N=50; 
% # time instants
T=1e6; 
% Is our data spikes (1) or flueresence traces (0) ?
isSpikes=0;

%% Generate network connectivity

con_pattern = {
'NIPS_network'
'combi'
'prob'        
'realistic'        
'realistic+1'        
'balanced'   
'balanced2'   
'common_input'   
'circular'
'rand'
'cluster'
'ringAtt'
'ringAtt2'
'block'};

% connectivity parameters (for details see "SetParams.m")
conn_type=con_pattern{4}; spar=0.1;inhib_frac=0.2;weight_dist='lognormal';seed_weights=1;weight_scale=1;N_stim=0;
% connectivity matrix
[W,~]=GetWeights(N,conn_type,spar,inhib_frac,weight_dist,seed_weights,weight_scale,N_stim,[],[]);
%bias
b=-3+0.1*randn(N,1); 


%% Generate spikes
addpath('GenerateSpikes');
% spiking parameters (for details see "SetParams.m")
T0=1e2;seed_spikes=1;neuron_type='logistic_with_history';N_stim=0; stim_type='pulses';timescale=1;s0=[],verbose=1;
% Get Spikes
S=GetSpikes(W,b,T,T0,seed_spikes,neuron_type,N_stim,stim_type,timescale,s0,verbose);

%% Figures 
figure(1)
subplot(2,2,1)
ma=max(W(:)); mi=min(W(:));
imagesc(W,[mi ma]); colorbar; title('W');  ylabel('#neuron') ; xlabel('#neuron')
title('Connectivity Pattern')

s = full(S);
subplot(2,2,2)
plot(s(1,1:100))
title('Example of neuron spiking signal (time bin 10 ms)')


subplot(2,2,3)
imagesc(s)
colorbar()
title('All the neurons')


corrmat = corrcoef(s');
corrmat(1:N+1:N*N) = 0;
subplot(2,2,4)
imagesc(corrmat)
title('Correlation Matrix')


figure(2) % figure2 is all about statistics
subplot(2,2,1)
ws = reshape(W,N*N,[]);
hist(ws, 100)
xlabel('weights')
ylabel('counts')

subplot(2,2,2)
cs = reshape(corrmat, [] ,1);
hist(cs, 100)
xlabel('correlations')
ylabel('counts')

% mean firing rates distribution
subplot(2,2,3)
fr = mean(s,2)/0.01;  % 10 ms 
fr_min = min(fr); fr_max = max(fr);
xbins = linspace(fr_min, fr_max, 20);
pp = hist(fr, xbins);
plot(xbins, pp, 'o--')
xlabel('firing rate')
ylabel('counts')

% population firing rates
subplot(2,2,4)
fr_pop = sum(s,1)/(N*0.01);


% shuffle each neuron
s_new = zeros(size(s));
for i=1:1:N
    ord = randperm(size(s,2));
    b = s(i,ord);
    s_new(i,:) = b;
end

fr_pop_rd = sum(s_new,1)/(N*0.01);

fr_pmin = min(min(fr_pop_rd), min(fr_pop)); 
fr_pmax = max(max(fr_pop_rd), max(fr_pop));

xbins = linspace(fr_pmin, fr_pmax, 10);
pfr_d = hist(fr_pop, xbins);
pfr_d_rd = hist(fr_pop_rd, xbins);

plot(xbins, pfr_d, 'r*-' , xbins, pfr_d_rd, 'bo-')

xlabel('Population firing rate')
ylabel('count')
legend('original','shuffle')
% save the data matrix
save('/Users/yanj11/data/data_logistic_spiking.mat', 's')
% save the weigtht matrix 
save('/Users/yanj11/data/data_logistic_spiking_w.mat', 'W')
figure(3)
imagesc(s_new)




