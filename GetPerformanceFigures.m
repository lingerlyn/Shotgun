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
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.1], [0.1 0.02], [0.1 0.02]);
T=2e6;
N=1e3;

% observations_ratios= [1,0.5,0.2,0.1,0.04,0.02];

K=5; %width of subplots
x_ticks={'R','C','Z','S'};
yticks=3; %number of y tick marks
ydigits=3; %number of  y digits
fontsize=12;
fontsize2=1.3*fontsize; 
if N==50
    dot_size=100; %for scatter plots
    observations_ratios= [1,0.2,0.1,0.04];
else
    dot_size=0.01;
    observations_ratios= [1,0.2,0.1];
end
    
L=length(observations_ratios);
a=length(observations_ratios);
b=3;

%%
fh=figure(1002);
set(fh,'units','normalized','outerposition',[0.2875    0.0521    0.7125    0.9479])
L_O=length(observations_ratios);
for ii=1:L_O
    load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(ii)) '_T=' num2str(T) '_Cavity_2.mat'],'W','EW','quality','centers','error_rates');
    
    %% Plot quality during convergence
    subplot(a,b,b*(ii-1)+1)
    iter=(1:length(quality(:,1)))/1000;
    [AX,H1,H2]=plotyy(iter,quality(:,2),iter,quality(:,3));
    if ii==L_O
        xlabel('1000 Iterations','fontsize',fontsize)
        legend('C','Z','location','SouthEast')
    end
%     set(get(AX(1),'Ylabel'),'String','C');
%     set(get(AX(2),'Ylabel'),'String','Z');

   

    for kk=1:2
       xlim(AX(kk),[-1e-1 max(iter)]);   
       ylim(AX(kk),[min(quality(:,kk+1)) 1]);    
       Lims=get(AX(kk),'Ylim');
       set(AX(kk),'YTick',ceil(linspace(Lims(1),Lims(2),yticks)*10^ydigits)/10^ydigits);

    end

    set(H2,'linestyle','--')
    label=[num2str(observations_ratios(ii))];
    if ii==1
        label={'$p_{\mathrm{obs}}$';'  ';'  '; label};
    end
    text(-0.4,0.5,label, 'Units', 'normalized','fontsize',fontsize2);    
    %% Plot quality as function of distance
    subplot(a,b,b*(ii-1)+2)
    [ quality_d,d_bins] = GetQualityDistance( W,EW,centers );
    [AX,H1,H2]=plotyy(d_bins,quality_d(:,1),d_bins,quality_d(:,3));
%     set(get(AX(1),'Ylabel'),'String','C');
%     set(get(AX(2),'Ylabel'),'String','Z');
    set(H2,'linestyle','--')
    if ii==L_O
        xlabel('Distance [$\mu$m]','fontsize',fontsize)
        legend('C','Z','location','NorthEast')
    end
    
    for kk=1:2
       xlim(AX(kk),[0 max(d_bins)])  
       ylim(AX(kk),[min(quality_d(:,2*(kk-1)+1)) 1]);    
       Lims=get(AX(kk),'Ylim');
       set(AX(kk),'YTick',ceil(linspace(Lims(1),Lims(2),yticks)*10^ydigits)/10^ydigits);
    end
    
    %% Plot ROCs
    [AUROC, ROC] = GetROCforW(error_rates);
    subplot(a,b,b*(ii-1)+3)
    [~,~,~,~,TPR_p,FPR_p,TPR_n,FPR_n] = GetWeightsErrors( W,EW );
    plot(squeeze(ROC(1,1,:)),squeeze(ROC(1,2,:)),'-b',squeeze(ROC(2,1,:)),squeeze(ROC(2,2,:)),'-r'...
        ,FPR_p,TPR_p,'xb',FPR_n,TPR_n,'xr','markers',12,'linewidth',2.5);
    if ii==L_O
        xlabel('FPR','fontsize',fontsize)
    end
    ylabel('TPR','fontsize',fontsize)
    ylim([0 1])
%     if ii==1
        legend(['E=' num2str(AUROC(1),2)],['I=' num2str(AUROC(2),2)],'location','southeast')
%     end    
    

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
    