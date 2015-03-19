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
subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.12], [0.1 0.05], [0.12 0.02]);
T=2e6;
N=50;

% observations_ratios= [1,0.5,0.2,0.1,0.04,0.02];

K=5; %width of subplots
x_ticks={'R','C','Z','S'};
fontsize=10;
fontsize2=1.3*fontsize; 
if N==50
    dot_size=100; %for scatter plots
    observations_ratios= [1,0.2,0.1,0.04];
else
    dot_size=0.01;
    observations_ratios= [1,0.2,0.1];
end
    
figure(1)
L=length(observations_ratios);
a=length(observations_ratios);
b=3;

%%
for ii=1:length(observations_ratios)
    load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(ii)) '_T=' num2str(T) '_Cavity.mat'],'W','EW','quality','centers','error_rates');
    
    %% Plot quality during convergence
    subplot(a,b,b*(ii-1)+1)
    iter=1:length(quality(:,1));
    [AX,H1,H2]=plotyy(iter,quality(:,2),iter,quality(:,3));
    xlabel('Iteration')
    set(get(AX(1),'Ylabel'),'String','C');
    set(get(AX(2),'Ylabel'),'String','Z');
    set(H2,'linestyle',':')
    xlim(AX(1),[0 max(iter)]);
    xlim(AX(2),[0 max(iter)]);
    label=[num2str(observations_ratios(ii))];
    if ii==1
        label={'$p_{\mathrm{obs}}$';'  ';'  '; label};
    end
    text(-0.5,0.5,label, 'Units', 'normalized','fontsize',fontsize2);    

    %% Plot quality as function of distance
    [AUROC, ROC] = GetROCforW(error_rates);
    subplot(a,b,b*(ii-1)+3)
    [~,~,~,~,TPR_p,FPR_p,TPR_n,FPR_n] = GetWeightsErrors( W,EW );
    plot(squeeze(ROC(1,1,:)),squeeze(ROC(1,2,:)),'-b',squeeze(ROC(2,1,:)),squeeze(ROC(2,2,:)),'-r'...
        ,FPR_p,TPR_p,'xb',FPR_n,TPR_n,'xr','markers',12);
    xlabel('FPR')
    ylabel('TPR')
    ylim([0 1])
    ylim([0 1])    
%     if ii==1
        legend(['E=' num2str(AUROC(1),2)],['I=' num2str(AUROC(2),2)],'location','southeast')
%     end    
    
    %% Plot quality as function of distance
    subplot(a,b,b*(ii-1)+2)
    [ quality_d,d_bins] = GetQualityDistance( W,EW,centers );
    [AX,H1,H2]=plotyy(d_bins,quality_d(:,1),d_bins,quality_d(:,3));
    set(get(AX(1),'Ylabel'),'String','C');
    set(get(AX(2),'Ylabel'),'String','Z');
    set(H2,'linestyle',':')
    xlabel('Distance [$\mu$m]')
    xlim(AX(1),[0 max(d_bins)]);
    xlim(AX(2),[0 max(d_bins)]);    
    ylim(AX(1),[min(quality_d(:,1)) 1]);
    ylim(AX(2),[min(quality_d(:,3)) 1]);    
end

target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript\Revision'
Export2Folder(['Performance_N=' num2str(N) '.eps'],target_folder) 


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
    