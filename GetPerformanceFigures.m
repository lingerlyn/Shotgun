clear all
% close all
clc

% Generate all figures for the paper
addpath('Results')
addpath('Misc')
addpath('GenerateSpikes');
addpath('Plotting')
SetDefaultGraphicSettings(0)
subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.08], [0.08 0.04], [0.12 0.02]);
T=2e6;
N=1e3;

% observations_ratios= [1,0.5,0.2,0.1,0.04,0.02];

K=5; %width of subplots
x_ticks={'R','C','Z','S'};
yticks=3; %number of y tick marks
ydigits=2; %number of  y digits
title_pos=[-0.2 0.98];

if N==50
    observations_ratios= [1,0.2,0.1,0.04];
else
    observations_ratios= [1,0.2,0.1];
end
    
L=length(observations_ratios);
a=length(observations_ratios);
b=3;
dark_green=[0 0.6 0.6];
dark_magenta=[0.6 0 0.6];

%%
fh=figure(1002);
% set(fh,'units','normalized','outerposition',[  0.5698    0.1204    0.4203    0.6574])
    units = 'centimeters';
    set(gcf, 'PaperUnits', units,'Units', units)           
    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'PaperPosition',[-0.0219   -0.2626   15   16])
    set(gcf, 'Position',[-0.0219   -0.2626   15   16])
L_O=length(observations_ratios);
for ii=L_O:-1:1
    load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(ii)) '_T=' num2str(T) '_Cavity_2.mat'],'W','EW','quality','centers','error_rates');
    
    %% Plot quality during convergence
    subplot(a,b,b*(ii-1)+1)
    iter=(1:length(quality(:,1)))/1000;
%     [AX,H1,H2]=plotyy(iter,quality(:,2),iter,quality(:,3));
    hp=plot(iter,quality(:,2),iter,quality(:,3));
    letter=['(' char(b*(ii-1)+'A') ')'] ;
    title(letter,'color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontweight','bold')
      if ii==L_O
        xlabel('1000 iterations')
      end
          if ii==1
        legend('C','Z','location','SouthEast')
    end
%     set(get(AX(1),'Ylabel'),'String','C');
%     set(get(AX(2),'Ylabel'),'String','Z');

   xlim([-1e-1 max(iter)]);   
   ylim([0 1]);    
   set(gca,'YTick',[0 0.25 0.5 0.75 1]);
%     for kk=1:2
%        xlim(AX(kk),[-1e-1 max(iter)]);   
% %        ylim(AX(kk),[min(quality(:,kk+1)) 1]);    
%         ylim(AX(kk),[0 1]);    
%        Lims=get(AX(kk),'Ylim');
% %        set(AX(kk),'YTick',ceil(linspace(Lims(1),Lims(2),yticks)*10^ydigits)/10^ydigits);
%         
%        set(AX(kk),'YTick',[0 0.2 0.4 0.6 0.8 1]);
% 
%     end
    
    set(hp(1),'color','k','linestyle','-')
    set(hp(2),'color',dark_green,'linestyle','--')
%     set(AX,{'ycolor'},{'k';dark_green});
    set(get(gca,'Ylabel'),'String',['p_{\rm{obs}}= ', num2str(observations_ratios(ii))]);

%% plot magnitude dependence
    subplot(a,b,b*(ii-1)+2)
    [ quality_m,m_bins] = GetQualityMagnitudeHist( W,EW);
%     [AX,H1,H2]=plotyy(m_bins,quality_m(:,1),m_bins,quality_m(:,2));
     hp=plot(m_bins,quality_m(:,1),m_bins,quality_m(:,2));
%     set(get(AX(1),'Ylabel'),'String','C');
%     set(get(AX(2),'Ylabel'),'String','Z');
    set(hp(2),'color',dark_magenta,'linestyle','-')
    set(hp(1),'color',dark_green,'linestyle','--')
%      set(AX,{'ycolor'},{dark_green;dark_magenta});
%     
    if ii==3
        xlabel('weights')       
    end
    if ii==1        
        lgnd=legend('Z','S','location','southeast');
    end
   xlim([-2.5 2.5]);   
   ylim([0 1]);    
   set(gca,'YTick',[0 0.25 0.5 0.75 1],'xtick',[-2 -1  0 1 2]');
%     for kk=1:2
%        xlim(AX(kk),[-1 1])  
%        ylim(AX(kk),[0 1]);    
%        Lims=get(AX(kk),'Ylim');
%        set(AX(kk),'YTick',[0 0.2 0.4 0.6 0.8 1],'xtick',[-1 -0.5 0 0.5 1]');
%     end
    letter=['(' char(b*(ii-1)+1+'A') ')'] ;
    title(letter,'color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontweight','bold')
    
       
    %% Plot ROCs
    [AUROC, ROC] = GetROCforW(error_rates);
    subplot(a,b,b*(ii-1)+3)
    [~,~,~,~,TPR_p,FPR_p,TPR_n,FPR_n] = GetWeightsErrors( W,EW );
    plot(squeeze(ROC(1,1,:)),squeeze(ROC(1,2,:)),'-.b',squeeze(ROC(2,1,:)),squeeze(ROC(2,2,:)),'-r'...
        ,FPR_p,TPR_p,'xb',FPR_n,TPR_n,'xr','markers',12);
    if ii==L_O
        xlabel('FPR')
    end
    ylabel('TPR')
    ylim([0 1])
%     if ii==1
        legend(['E=' num2str(AUROC(1),2)],['I=' num2str(AUROC(2),2)],'location','southeast')
%     end    

    letter=['(' char(b*(ii-1)+2+'A') ')'] ;
    title(letter,'color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontweight','bold')

end

target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript\Revision2'
Export2Folder(['Fig7.eps'],target_folder) 

 %% Plot quality as function of distance
% fh=figure(1003);
% subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.1], [0.15 0.05], [0.2 0.15]);
% 
% set(fh,'units','normalized','outerposition',[ 0.3245    0.4380    0.1875    0.4667])
% 
% for ii=1:L_O
%     load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(ii)) '_T=' num2str(T) '_Cavity_2.mat'],'W','EW','quality','centers','error_rates');
%      
%     subplot(3,1,ii)
%     [ quality_d,d_bins] = GetQualityDistance( W,EW,centers );
%     [AX,H1,H2]=plotyy(d_bins,quality_d(:,1),d_bins,quality_d(:,3));
% %     set(get(AX(1),'Ylabel'),'String','C');
% %     set(get(AX(2),'Ylabel'),'String','Z');
%     set(H2,'linestyle','--')
%     
%     if ii==3
%         xlabel('Distance [$\mu$m]','fontsize',fontsize)
%     end
%     if ii==2
%             lgnd=legend('C','Z','location','NorthEast')
%         set(lgnd,'fontsize',fontsize)
%     end
%     for kk=1:2
%        xlim(AX(kk),[0 max(d_bins)])  
%        ylim(AX(kk),[min(quality_d(:,2*(kk-1)+1)) 1]);    
%        Lims=get(AX(kk),'Ylim');
%        set(AX(kk),'YTick',ceil(linspace(Lims(1),Lims(2),yticks)*10^ydigits)/10^ydigits,'xlim',[0 150]);
%     end
% 
%      set(get(AX(1),'Ylabel'),'String',['$p_{\mathrm{obs}}$=', num2str(observations_ratios(ii))],'fontsize',fontsize+1);
% 
% end
% 
% target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript\Revision2'
% Export2Folder(['Distance_Dependence.jpg'],target_folder) 
%     