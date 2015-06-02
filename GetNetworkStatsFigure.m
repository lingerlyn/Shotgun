% clear all
% close all
clc

% Generate all figures for the paper
addpath('Results')
addpath('Misc')
addpath('Plotting')
addpath('GenerateSpikes');
addpath('GenerateConnectivity')
SetDefaultGraphicSettings(0)

fontsize=12;
linewidth=2;
title_pos=[-0.35, 1.04];
title_pos2=[-0.06, 1.005];
 set(0,'defaultlinelinewidth',2)

set(0,'DefaultTextInterpreter', 'tex');
% set(0,'DefaultAxesFontSize',fontsize-2)
a=3; b=4;
subplot = @(p) subtightplot (a,b, p, [0.06 0.08], [0.06 0.05], [0.1 0.02]);

params=SetParams();
if params.spike_gen.T>1e4
    params.spike_gen.T=1e4;
end
T_view=3e2;
dt=1e-2; %100 Hz imaging frame rate
regenerate_data=1;

h=figure(1001);
% set(h,'units','normalized','outerposition',[ 0.5240    0.0528    0.4714    0.8574])
units = 'centimeters';
set(gcf, 'PaperUnits', units,'Units', units)           
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperPosition',[-0.0219   -0.2626   17   18])
set(gcf, 'Position',[-0.0219   -0.2626   17   18])
%%  W

if regenerate_data

    [N,spar,inhib_frac,weight_dist,bias,seed_weights, weight_scale, conn_type,N_stim,target_rates]=v2struct(params.connectivity);
    [W,~]=GetWeights(N,conn_type,spar,inhib_frac,weight_dist,seed_weights,weight_scale,N_stim,params.spike_gen.stim_type,params.sbm);
    if ~isempty(target_rates)
        bias=SetBiases(W,target_rates,params.spike_gen);
    end
    % bias=bias-sqrt(-bias/2).*randn(N,1);
    mi=min(W(:));
    ma=max(W(:));
    if mi>-ma
        mi=-ma;
    elseif ma<-mi
        ma=-mi;
    end
end

%%
% subplot([9 10 11 13 14  15])
% blue_red_map2= darkb2r(mi,ma);
% imagesc(W,[mi ma]);
% xlabel('Cell index');
% ylabel('Cell index');
% colormap(blue_red_map2)
% freezeColors
% cbar=colorbar('location','eastoutside');
% cpos = get(cbar,'position');
% cpos = cpos + [0.07,0,0,0]; 
% set(cbar,'position',cpos)
% title('(B)','color', 'k', 'fontweight', 'bold', 'Units', 'normalized','position',title_pos3,'fontsize',fontsize)
% cbfreeze(cbar)
%%
subplot(4)
[hist_W,bins]=hist(W(~~W(:)),100);
hist_W=hist_W/sum(hist_W);
plot(bins,hist_W, 'k')
xlabel('weights')
ylabel('count');
title('(B)','color', 'k', 'fontweight', 'bold', 'Units', 'normalized','position',title_pos,'fontsize',fontsize)
xlim([mi ma])
box('off')
%% Plot Correlations
load('Run_N=1000_obs=1_T=2000000_Cavity.mat','Cxy','Cxx');
subplot(8)
vars=sqrt(diag(Cxx));
corr_delay=Cxy./(vars*vars');
corr_delay(eye(size(corr_delay))>0.5)=0; 
% corr=Cxx./(vars*vars');
% imagesc(corr_delay)
[hist_values,bins]=hist(corr_delay(~~(corr_delay(:))),30);
% C=sum(hist_values);
% hist_values=hist_values/C;
semilogy(bins,hist_values, 'k')
xlim([-0.1 bins(end)])
% ylim([1e-6 1])
% set(gca,'ytick',[1e-6 1e-4 1e-2 1])
ylim([1 1e6])
set(gca,'ytick',[1 1e2 1e4 1e6])
ylabel('count');
xlabel('correlation');
title('(C)','color', 'k', 'fontweight', 'bold', 'Units', 'normalized','position',title_pos,'fontsize',fontsize)
large_correlations=sum(hist_values(bins>=0.1))/sum(hist_values)
negative_correlations=sum(hist_values(bins<0))/sum(hist_values)
box('off')
%% Spikes - generate data
if regenerate_data
    addpath('GenerateSpikes');
    [T,T0,sample_ratio,sample_type,seed_spikes,seed_sample,N_stim,stim_type, neuron_type,timescale,obs_duration]=v2struct(params.spike_gen);
    s0=[];verbose=1;
    spikes=GetSpikes(W,bias,T,T0,seed_spikes,neuron_type,N_stim,stim_type,timescale,s0,verbose);
    for ii=1:N
        spikes_shuffled(ii,:)=spikes(ii,randperm(T)); %#ok
    end
    rho=corr(spikes);
end
%% Spikes - plot
subplot([1 2 3   5 6 7  9 10 11 ])

tt=(1:T_view)*dt; nn=1:N;
imagesc(tt,nn,1-spikes(nn,1:T_view))
% hbar=colorbar
colormap('gray')
xlabel('time [sec]');
ylabel('neuron');
line([-1 2*tt(end)],1000*[1 1],'color','red','linewidth',2)
freezeColors
% cbfreeze(hbar)
title('(A)','color', 'k', 'fontweight', 'bold', 'Units', 'normalized','position',title_pos2,'fontsize',fontsize)
% "Burst size" distribution
box('off')
%%
subplot(12)
L=20;
[hist_spikes,bins]=hist(mean(spikes,1),L);
[hist_shuffled,bins_shuffled]=hist(mean(spikes_shuffled,1),L);
semilogy(bins/dt,hist_spikes,'ob',bins_shuffled/dt,hist_shuffled,'*r','MarkerSize',4)
ylabel('count');
xlabel('firing rate [Hz]');
% h=legend('Data','Shuffled');
% set(h,'fontsize',10);
title('(D)','color', 'k', 'fontweight', 'bold', 'Units', 'normalized','position',title_pos,'fontsize',fontsize)
box('off')
% figure
% rho_vec=rho(eye(N)<0.5);
% hist(rho_vec,100);

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
% set(gca, ...
%   'Box'         , 'off'     , ...
%   'TickDir'     , 'out'     , ...
%   'TickLength'  , [.05 .05] , ...
%   'XMinorTick'  , 'on'      , ...
%   'YMinorTick'  , 'on'      , ...
%   'YGrid'       , 'on'    )  

%%

target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript\Revision2';
Export2Folder(['Fig2.eps'],target_folder) 