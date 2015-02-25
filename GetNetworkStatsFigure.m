clear all
% close all
clc

% Generate all figures for the paper
addpath('Results')
addpath('Misc')
addpath('GenerateSpikes');
addpath('GenerateConnectivity')


set(0,'DefaultTextInterpreter', 'latex');
set(0,'DefaultAxesFontSize',13)
a=1; b=4;
subplot = @(p) subtightplot (a,b, p, [0.06 0.08], [0.1 0.06], [0.05 0.02]);

params=SetParams();
if params.spike_gen.T>1e4
    params.spike_gen.T=1e4;
end
T_view=1e3;
dt=1e-2; %100 Hz imaging frame rate

h=figure(1001);
set(h,'units','normalized','outerposition',[0 0 1.2 0.6])
%%  W

[N,spar,inhib_frac,weight_dist,bias,seed_weights, weight_scale, conn_type,N_stim,target_rates]=v2struct(params.connectivity);
[W,~]=GetWeights(N,conn_type,spar,inhib_frac,weight_dist,seed_weights,weight_scale,N_stim,params.spike_gen.stim_type,params.sbm);
if ~isempty(target_rates)
    bias=SetBiases(W,target_rates,params.spike_gen);
end
% bias=bias-sqrt(-bias/2).*randn(N,1);

% subplot(7)
% imagesc(W);
% xlabel('cell number');
% ylabel('cell number');
% % colormap('jet')
% freezeColors
% h=colorbar;
% cbfreeze(h)
% 
% subplot(8)
% [hist_W,bins]=hist(W(~~W(:)),1000);
% plot(bins,hist_W)
% xlabel('weight');
% ylabel('count');

%% Spikes - generate data
addpath('GenerateSpikes');
[T,T0,sample_ratio,sample_type,seed_spikes,seed_sample,N_stim,stim_type, neuron_type,timescale,obs_duration]=v2struct(params.spike_gen);
s0=[];verbose=1;
spikes=GetSpikes(W,bias,T,T0,seed_spikes,neuron_type,N_stim,stim_type,timescale,s0,verbose);
for ii=1:N
    spikes_shuffled(ii,:)=spikes(ii,randperm(T)); %#ok
end
rho=corr(spikes);
%% Spikes - plot
subplot([1 2 3])

tt=(1:T_view)*dt; nn=1:N;
imagesc(tt,nn,spikes(:,1:T_view))
colorbar
colormap('gray')
xlabel('Time [sec]');
ylabel('Cell index');
line([-1 2*tt(end)],1000*[1 1],'color','red','linewidth',2)
% "Burst size" distribution

subplot(4)
L=20;
[hist_spikes,bins]=hist(mean(spikes,1),L);
[hist_shuffled,bins_shuffled]=hist(mean(spikes_shuffled,1),L);
semilogy(bins/dt,hist_spikes,'ob',bins_shuffled/dt,hist_shuffled,'*r')
ylabel('Frequency');
xlabel('Firing rate [Hz]');
legend('Data','Shuffled');

figure
rho_vec=rho(eye(N)<0.5);
hist(rho_vec,100);

%% Observations
% addpath('GenerateSpikes');
% t_start=0; 
% observations=SampleSpikes(N,T,sample_ratio,sample_type,obs_duration,N_stim,seed_sample,t_start);
% 
% subplot([4 5 6])
% imagesc(observations(:,1:T_view));
% xlabel('time');
% ylabel('cell number');
% 
% XXn=observations*observations'/T;
% subplot(3+b)
% imagesc(XXn);
% ylabel('cell number');
% ylabel('cell number');

%% Autocovariance
% c_all=0;
% for ii=1:N
%     [c,lags]=xcov(full(spikes(ii,:)));
%     c_all=c_all+c;
% end
% plot(lags,c_all/N);

%% firing rate distributions
% W_temp=W;
% N=size(W,1);
% W_temp(eye(N)>0.5)=0;
% ind_excitatory=(sum(W_temp,1)>=0)';
% ind_inhibitory=(sum(W_temp,1)<0)';
% 
% spikes_exc=spikes(ind_excitatory,:);
% spikes_inh=spikes(ind_inhibitory,:);
% 
% L=30;
% bins=logspace(-4,0,L)
% figure
% [hist_exc,bins]=hist(mean(spikes_exc,2),bins);
% semilogx(bins,hist_exc)
% figure
% [hist_inh,bins]=hist(mean(spikes_inh,2),bins);
% semilogx(bins,hist_inh)

%%
target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript\Revision';
Export2Folder(['NetworkStats.eps'],target_folder) 