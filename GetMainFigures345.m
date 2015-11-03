clear all
% close all
clc


N=50; %Number of neurons. Set N=50 for figures 3&4, N=1e3 for figure 5
isLIF=0; %are we using LIF? %Set isLIF=0 for figures 3&5, isLIF=1 for figure 4

% Generate all figures for the paper
addpath('Results')
addpath('Misc')
addpath('Plotting')
addpath('GenerateSpikes');
% set(0,'DefaultTextInterpreter', 'tex');
% set(0,'DefaultAxesFontSize',10)
SetDefaultGraphicSettings(0)
T=2e6;
T_view=3e2;
dt=1e-2; %100 Hz imaging frame rate


figure(1098)
if N==50
    subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.08], [0.06 0.05], [0.04 0.01]);
    observations_ratios= [1,0.2,0.1,0.04];
    dot_size=3.8;
    ind_array=1:N;
    shiftbar=[0.075,0,0,0];
    shiftpos= [0.025,0,0,0]; %shift left column becauase of (aspect ratio)
    units = 'centimeters';
    set(gcf, 'PaperUnits', units,'Units', units)           
    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'PaperPosition',[-0.0219   -0.2626   17   18])
    set(gcf, 'Position',[-0.0219   -0.2626   17   18])
elseif N==1e3;
    subplot = @(m,n,p) subtightplot (m, n, p, [0.008 0.08], [0.08 0.03], [0.06 0.005]);
    observations_ratios= [1,0.2,0.1];
    dot_size=3.8;
    M=50;
    ind_array=randi(N,M,1);
    shiftbar=[0.1,0,0,0];
    shiftpos= [0.015,0,0,0]; %shift left column becauase of (aspect ratio)
    units = 'centimeters';
    set(gcf, 'PaperUnits', units,'Units', units)           
    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'PaperPosition',[-0.0219   -0.2626   17   16])
    set(gcf, 'Position',[-0.0219   -0.2626   17   16])
end

fontsize2=12; 
L_bins=30; %number of bins in the weight histogram
title_pos=[-0.25 0.85]; % title position for colmun 1
title_pos2=[0.1 0.85]; % title position for colmun 2
title_pos3=[0.1 0.85]; % title position for colmun 3
title_pos4=[-0.7 0.85]; % title position for colmun 4
%% Figure 1 - Toy model 
% In this figure:
% N=50, 
% T=2e5 (~1 hour experiment, assuming 200ms bins)
% firing rate~0.2   (~1Hz firing rate)
% balanced network (no Dale's law yet)

L=length(observations_ratios);
K=7; %width of subplots
b=4; %number of columns
x_ticks={'R','C','Z','S'};

blue_red_map2= darkb2r(-1,1);

%% Plot True weight matrix

% set(gcf,'outerposition',[961   ,       49 ,        960    ,    1032])


if isLIF==1
    LIF_str='_LIF';
else
    LIF_str=[];
end
load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(2)) '_T=' num2str(T) LIF_str '_Cavity.mat'],'W','W_full','bias_full','EW','bias','rates','params','centers','quality');
W_temp=W(ind_array,ind_array);
mi=min(W_temp(:));
ma=max(W_temp(:));
if mi>-ma
    mi=-ma;
elseif ma<-mi
    ma=-mi;
end
% mi=-1;ma=1;
mi2=min(W(:))*1.2;
ma2=max(W(:))*1.2;

subplot(L+1,K,[1 2])    
imagesc(W(ind_array,ind_array),[mi ma]);
set(gca,'DataAspectRatio',[1 1 1],'xtick',[])
ylabel('true W','fontsize',fontsize2)
colormap(blue_red_map2)
freezeColors
originalSize = get(gca, 'Position');
h=colorbar('location','eastoutside');
h=cbfreeze(h);
cpos = get(h,'position');
cpos = cpos + shiftbar; 
set(h,'position',cpos,'ylim', [mi ma])
set(gca, 'Position',originalSize + shiftpos)
title('(A)','color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontweight','bold')



      

% subplot(L+1,K,[3 4])
% [hist_W,bins]=hist(W(~~W(:)),L_bins);
% plot(bins,hist_W)
% % xlabel('W');
% ylabel('count');

% ii=0;
%  subplot(L+1,K,K*ii+[5 6])
% A_ind=linspace(mi2,ma2,100);
% plot(A_ind,A_ind,'r-');
% hold all
% W_temp=W(ind_array,:);
% plot(W_temp(:),W_temp(:),'b.','MarkerSize',dot_size)
% %     cloudPlot(W(:),EW(:),[mi ma mi ma],1,[30 30])
% axis([mi2 ma2 mi2 ma2])
% hold off
% %     legend('x=y','EW','EW2')
% % xlabel('W','fontsize',fontsize2)
% ylabel('W','fontsize',fontsize2)

% subplot(L+1,K,K*ii+[3 4])
% mi3=-2;%min([W(:) ;EW(:)]);
% ma3=2;%max([W(:) ;EW(:)]);
% bins=linspace(mi3,ma3,L_bins);
% %     hist_W=hist(W(~~W(:)),bins);
% %     hist_EW=hist(EW(~~EW(:)),bins);
% hist_W=hist(W(:),bins);
% semilogy(bins,hist_W,'bo','linewidth',2);
% xlim([mi3 ma3]);
% xlabel('W','fontsize',fontsize2)
% ylabel('count','fontsize',fontsize2)
%     if ii==1
%         h=legend('W','$\hat{W}$')
%         set(h,'Interpreter','latex','fontsize',fontsize,'location','northwest','orientation','horizontal')
%     end

% subplot(L+1,K,K*ii+7)
% [R,correlation, zero_matching,sign_matching] = GetWeightsErrors( W,W );
% % AUROC=[1 1];
% % bar( [R,correlation, zero_matching,sign_matching,AUROC(1),AUROC(2)] );    
% bar( [R,correlation, zero_matching,sign_matching] );    
% ylim([0 1])
% ylabel('quality','fontsize',fontsize)
% set(gca, 'XTickLabel', x_ticks,'fontsize',fontsize);


%%
for ii=L:-1:1
    load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(ii)) '_T=' num2str(T) LIF_str '_Cavity.mat'],'W','EW','rates');
        %     high_firing_rate=rates>0.05;
%     W=W(high_firing_rate,high_firing_rate);
%     EW=EW(high_firing_rate,high_firing_rate);
       
    subplot(L+1,K,K*ii+[1 2])    
    imagesc(EW(ind_array,ind_array),[mi ma]);
    set(gca,'DataAspectRatio',[1 1 1],'xtick',[])
    if ii~=L
        set(gca,'xtick',[])
    else
        set(gca,'xtick',[ 10 20 30 40 50])
    end
    h=colorbar;
    set(h, 'ylim', [mi ma])
    ylabel('true W','fontsize',fontsize2)
    colormap(blue_red_map2)
    freezeColors
%     cbfreeze(h)
    colorbar('off');
    ylabel(['p_{\rm{obs}}= ' num2str(observations_ratios(ii))],'fontsize',fontsize2)
    colormap('jet')
    pos = get(gca,'position');
    pos = pos +  shiftpos; 
    set(gca, 'Position',pos)
    letter=['(' char(b*(ii-1)+1+'A') ')'] ;
    title(letter,'color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontweight','bold')

    

    subplot(L+1,K,K*ii+[5 6])
    A_ind=linspace(mi2,ma2,100);
    plot(A_ind,A_ind,'k--','linewidth',1);
    hold all
    EW_temp=EW(ind_array,:);
    W_temp=W(ind_array,:);
    plot(W_temp(:),EW_temp(:),'b.','MarkerSize',dot_size)
%     cloudPlot(W(:),EW(:),[mi ma mi ma],1,[30 30])
    axis([mi2 ma2 mi2 ma2])
    hold off
%     legend('x=y','EW','EW2')
if ii==L
    xlabel('true weights','fontsize',fontsize2)
else
    set(gca,'xtick',[]);
end
    ylabel('inferred weights','fontsize',fontsize2)
letter=['(' char(b*(ii-1)+3+'A') ')'] ;
    title(letter,'color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos3,'fontweight','bold')
    
    subplot(L+1,K,K*ii+[3 4])
    mi3=-2;%min([W(:) ;EW(:)]);
    ma3=2;%max([W(:) ;EW(:)]);
    bins=linspace(mi3,ma3,L_bins);
%     hist_W=hist(W(~~W(:)),bins);
%     hist_EW=hist(EW(~~EW(:)),bins);
    hist_W=hist(W(:),bins);
    hist_EW=hist(EW(:),bins);
    hs=semilogy(bins,hist_W,'bo',bins,hist_EW,'r.');
    for kk=1:2
        set(hs(kk),'MarkerSize',5)
    end
    xlim([mi3 ma3]);
    max_hist=log10(max(hist_W(:)));
    ylim([1 3*10^max_hist])
    set(gca,'ytick',10.^(0:floor(max_hist)));
    if ii==L
        xlabel('weights','fontsize',fontsize2)
    else
        set(gca,'xtick',[]);
    end
    ylabel('count','fontsize',fontsize2)
    letter=['(' char(b*(ii-1)+2+'A') ')'] ;
    title(letter,'color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos2,'fontweight','bold')
%     if ii==1
%         h=legend('W','$\hat{W}$')
%         set(h,'Interpreter','latex','fontsize',fontsize,'location','northwest','orientation','horizontal')
%     end
    
   subplot(L+1,K,K*ii+7)
   [R,correlation, zero_matching,sign_matching] = GetWeightsErrors( W,EW );
%    [AUROC, ROC] = GetROCforW(quality)
    % bar( [R,correlation, zero_matching,sign_matching,AUROC(1),AUROC(2)] );    
    bar( [R,correlation, zero_matching,sign_matching] ,'k');    
    ylim([0 1])
    ylabel('quality')
    set(gca,'ytick',[0.25 0.5 0.75 1]);
    if ii==L
        set( gca, 'XTickLabel', x_ticks); 
    else
        set( gca, 'XTickLabel',[]); 
    end
    letter=['(' char(b*(ii-1)+4+'A') ')'] ;
    title(letter,'color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos4,'fontweight','bold')


end
% 
target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript\FinalProduction';
if isLIF
    Export2Folder([ 'Fig4.tif'],target_folder) 
else
    if N==50
      Export2Folder([ 'Fig3.tif'],target_folder) 
    elseif N==1e3
      Export2Folder([ 'Fig5.tif'],target_folder) 
    end
end

% Export2Folder(['Sparsity' LIF_str '_N=' num2str(N) '.eps'],target_folder) 