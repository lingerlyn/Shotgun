clear all
close all
clc


% Generate all figures for the paper
addpath('Results')
addpath('Misc')
% addpath('GenerateSpikes');
% addpath('GenerateConnectivity')
addpath('Plotting')
SetDefaultGraphicSettings( )

% set(0,'DefaultTextInterpreter', 'latex');
% set(0,'DefaultAxesFontSize',13)

fontsize=12;
x_ticks={'R','C','Z','S'};
% set(0,'defaultlinelinewidth',1)

dt=1e-2; %100 Hz imaging frame rate

N=50;
T=5e5;
observations_ratios=[1 0.1];
title_pos=[-0.15 1.05]
title_pos2=[-0.3  1.05];

%% Matrix reconstruction
h=figure(1002);
% set(h,'units','normalized','outerposition',[0 0.2 0.3 0.6])
    units = 'centimeters';
    set(gcf, 'PaperUnits', units,'Units', units)           
    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'PaperPosition',[-0.0219   -0.2626   12   12])
    set(gcf, 'Position',[-0.0219   -0.2626   12   12])

a=2; b=3;
subplot = @(p) subtightplot (a,b, p, [0.1 0.12], [0.1 0.1], [0.15 0.02]);
for ii=1:length(observations_ratios)
    load(['Run_N=' num2str(N) '_obs=' num2str(observations_ratios(ii)) '_T=' num2str(T) '_EM.mat'],'W','EW','EW2','EW3','params');
    
        subplot([1 2]+(ii-1)*b)
        mi=min([W(:); EW(:)])*1.2;
        ma=max([W(:); EW(:)])*1.2;
        A_ind=linspace(mi,ma,100);
        hp=plot(A_ind,A_ind,'k--','linewidth',1);
        hold all
        plot(W(:),EW(:),'.b');
        hold all
        if ii==1
            plot(W(:),EW2(:),'.r');
        else
            plot(W(:),EW3(:),'.m');
        end
        hold all
        axis([mi ma mi ma])
        xlabel('true weights','fontsize',fontsize)
        ylabel({['p_{\rm{obs}}= ' num2str(observations_ratios(ii))] ,' inferred weights'},'fontsize',fontsize)
        letter=['(' char(2*(ii-1)+'A') ')'] ;
        title(letter,'color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontweight','bold')


subplot(3+(ii-1)*b)
[R,correlation, zero_matching,sign_matching] = GetWeightsErrors( W,EW ); 
[R2,correlation2, zero_matching2,sign_matching2] = GetWeightsErrors( W,EW2 ); 
[R3,correlation3, zero_matching3,sign_matching3] = GetWeightsErrors( W,EW3 ); 
if ii==1
    h=bar(1:4,[ [R,correlation, zero_matching,sign_matching];  [R2,correlation2, zero_matching2,sign_matching2]]',2);  
    set(h(1),'FaceColor','b');
    set(h(2),'FaceColor','r');
else
    h=bar(1:4,[ [R,correlation, zero_matching,sign_matching];  [R2,correlation2, zero_matching2,sign_matching2] ; [R3,correlation3, zero_matching3,sign_matching3] ]',2); 
    set(h(1),'FaceColor','b');
    set(h(2),'FaceColor','c');
    set(h(3),'FaceColor','m');
end

ylim([0 1])
ylabel('quality','fontsize',fontsize)
set(gca, 'XTickLabel', x_ticks,'fontsize',fontsize);
letter=['(' char(2*(ii-1)+1+'A') ')'] ;
title(letter,'color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos2,'fontweight','bold')
end
target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript\FinalProduction';
Export2Folder(['Fig9.tif'],target_folder) 