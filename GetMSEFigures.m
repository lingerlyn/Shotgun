clear all
% close all
clc

% Generate all figures for the paper
addpath('Results')
addpath('Misc')
addpath('GenerateSpikes');
set(0,'DefaultTextInterpreter', 'latex');
set(0,'DefaultAxesFontSize',11)
set(0,'defaultlinelinewidth',2)
subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.15], [0.1 0.05], [0.1 0.01]);
T=2e6;
N=50;

% observations_ratios= [1,0.5,0.2,0.1,0.04,0.02];
observations_ratios= [1,0.2,0.1,0.04];
L=length(observations_ratios);
K=5; %width of subplots
x_ticks={'R','C','Z','S'};
fontsize=10;
fontsize2=1.3*fontsize; 
if N==50
    dot_size=100; %for scatter plots
else
    dot_size=0.01;
end
    
figure(1)
a=length(observations_ratios);
b=2;

%%
for ii=1:length(observations_ratios)
    load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(ii)) '_T=' num2str(T) '_Cavity.mat'],'W','EW','quality','centers');
    
    %% Plot quality during convergence
    subplot(a,b,b*(ii-1)+1)
    plot(quality(:,1))
    xlabel('Iteration')
    ylabel('C (correlation)')
    xlim([0 1e4])
    text(5e3,0.2,['$p_{\mathrm{obs}}=$' num2str(observations_ratios(ii))],'fontsize',fontsize2)

    %% Plot quality as function of distance
    subplot(a,b,b*(ii-1)+2)
    [ quality_d,d_bins] = GetQualityDistance( W,EW,centers );
    plot(d_bins,quality_d(:,1))
    xlabel('Distance [$\mu$m]')
    ylabel('C (correlation)')
    xlim([0 60])
    mi=min(quality_d(:,1));
    ylim([mi 1])
    text(10,mi+0.2*(1-mi),['$p_{\mathrm{obs}}=$' num2str(observations_ratios(ii))],'fontsize',fontsize2)
end

target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript\Revision'
Export2Folder(['Convergence_figure.eps'],target_folder) 


%% Plot quality during convergence
% figure
% for kk=1:4
%     subplot(4,1,kk)
%     plot(quality(:,kk))
% end
% 
% %% Plot quality as function of distance
% figure
% for kk=1:4
%     [ quality_d,d_bins] = GetQualityDistance( W,EW,centers );
%     subplot(4,1,kk)
%     plot(d_bins,quality_d(:,kk))
% end
    