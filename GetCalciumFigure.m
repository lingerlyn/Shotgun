clear all
close all
clc


% Generate all figures for the paper
addpath('Results')
addpath('Misc')
% addpath('GenerateSpikes');
% addpath('GenerateConnectivity')

set(0,'DefaultTextInterpreter', 'latex');
set(0,'DefaultAxesFontSize',13)

fontsize=13;
x_ticks={'R','C','Z','S'};
set(0,'defaultlinelinewidth',1)

dt=1e-2; %100 Hz imaging frame rate

N=50;
T=2e6;
observations_ratios=[1 0.2 0.1];
title_pos=[-0.15 0.9];

load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(1)) '_T=' num2str(T) '_Cavity_CalciumObs.mat'],'W','EW', 'Spike_rec_correlation','sample_traces');


%% spike reconstructions
h=figure(1001);
set(h,'units','normalized','outerposition',[0 0.2 0.7 0.4])
a=1; b=4;
subplot = @(p) subtightplot (a,b, p, [0.06 0.1], [0.2 0.06], [0.05 0.02]);

% T_view=size(sample_traces.Y,2);
T_view=4e3;
tt=(1:T_view)*dt;
ii=6;
subplot([1 2 3])
plot(tt,sample_traces.Y(ii,1:T_view)/max(sample_traces.Y(ii,1:T_view))+0.2,'k')
xlim([0 max(tt)])
hold all
plot(tt,sample_traces.spikes(ii,1:T_view),'ob',tt,sample_traces.true_spikes(ii,1:T_view),'.r')
ylim([0.1 1.1]);
xlabel('Time [Sec]','fontsize',fontsize);
set(gca,'ytick',[]);
% title('A','color', 'k', 'fontweight', 'bold', 'Units','normalized', 'interpreter','none','position',[-0.03 0.95],'fontsize',fontsize)

% spike reconstruction correlations
T_view=size(sample_traces.Y,2);
N_view=size(sample_traces.Y,1);
tt=(1:T_view)*dt;
ii=6;
subplot(4)
plot(Spike_rec_correlation,'linewidth',2)

xlabel('Cell index','fontsize',fontsize)
ylabel('Correlation','fontsize',fontsize);
xlim([1 N_view]);
% title('B','color', 'k', 'fontweight', 'bold', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontsize',fontsize)

target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript\Revision';
Export2Folder(['Calcium_spikes.eps'],target_folder) 
%% Matrix reconstruction
h=figure(1002);
set(h,'units','normalized','outerposition',[0 0.2 0.3 0.8])
a=3; b=3;
subplot = @(p) subtightplot (a,b, p, [0.08 0.12], [0.1 0.06], [0.15 0.02]);
for ii=1:length(observations_ratios)
    load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(ii)) '_T=' num2str(T) '_Cavity_CalciumObs.mat'],'EW','W');
    EW_c=EW;
    W_c=W;

    [R_c,correlation_c, zero_matching_c,sign_matching_c] = GetWeightsErrors( W,EW_c ); 


load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(ii)) '_T=' num2str(T) '_Cavity_noCalcium.mat'],'W','EW');

subplot([1 2]+(ii-1)*b)
mi=min(W(:))*1.2;
ma=max(W(:))*1.2;
A_ind=linspace(mi,ma,100);
plot(A_ind,A_ind,'g-')
hold all
plot(W(:),EW(:),'.b');
hold all
plot(W_c(:),EW_c(:),'.r');
axis([mi ma mi ma])
xlabel('W','fontsize',fontsize)
ylabel({['$p_{\mathrm{obs}}=$' num2str(observations_ratios(ii))] ,'\, \quad$\hat{W}$'},'fontsize',fontsize)
% title('E','color', 'k', 'fontweight', 'bold', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontsize',fontsize)

subplot(3+(ii-1)*b)
[R,correlation, zero_matching,sign_matching] = GetWeightsErrors( W,EW ); 
bar(1:4,[ [R,correlation, zero_matching,sign_matching];  [R_c,correlation_c, zero_matching_c,sign_matching_c]]',2);    
ylim([0 1])
ylabel('Quality','fontsize',fontsize)
set(gca, 'XTickLabel', x_ticks,'fontsize',fontsize);
% title('F','color', 'k', 'fontweight', 'bold', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontsize',fontsize)
end
target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript\Revision';
Export2Folder(['Calcium_weights.eps'],target_folder) 