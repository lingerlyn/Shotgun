clear all
% close all
clc

%%

addpath('Results/ParametersScanResults')
addpath('Plotting')

SetDefaultGraphicSettings(1)
h=figure(1001);

% colormap('jet' )
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.11], [0.18 0.04], [0.12 0.01]);
labelfontsize=12;
t_start=4;
title_pos=[-0.18 0.98];
% set(h,'units','normalized','outerposition',[0 0.02 0.6 0.8])
units = 'centimeters';
set(gcf, 'PaperUnits', units,'Units', units)           
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperPosition',[4.5073   10.3273   16    14])
set(gcf, 'Position',[4.5073   10.3273   16   12])

for kk=1:4
    if kk==1
        load('scan_Pobs.mat');
        % cut out regions that did not converge
%         correlation_mat (:,1:2)=[];
        y_name='p_{\rm{obs}} - observed fraction';
        T_str=2.^((t_start-1):10);
         y_str =2.^(-(0:5));
         L_T=length(T_str);
         correlation_mat = correlation_mat(:,t_start:end);
    elseif kk==2
        load('scan_N.mat');
        
        y_str=2.^(11:-1:5);
        T_str=2.^((t_start-1):10);        
        ind_n=length(y_str):-1:1;
        L_T=length(T_str);
        correlation_mat = correlation_mat(ind_n,t_start:end);
        y_name='N - number of neurons';
    elseif kk==3
        load('scan_target_N=500.mat');
        y_str=[0.32 0.16 0.08 0.04 0.02 0.01]*100;
        T_str=2.^((t_start-1):10);
        L_T=length(T_str);
%         correlation_mat = correlation_mat(:,t_start:end,2);
        correlation_mat = correlation_mat(:,t_start-2:end);
        y_name='m - firing rate [Hz]';
    elseif kk==4
         load('scan_spar_N=500.mat');
        y_str=[0.4 0.2 0.1 0.05 0.025];
        T_str=2.^((t_start-1):10);
        L_T=length(T_str);
%         correlation_mat = mean(correlation_mat(:,t_start:end,1:4),3);
        correlation_mat = correlation_mat(:,t_start-2:end);
        y_name='p_{ 0} - connection sparsity';
    elseif kk==5
         load('scan_wscale.mat');
        y_str=[1 0.5 0.25 0.125];
        R_squared_mat(1,:)=[];
        correlation_mat(1,:)=[];
        ZM_mat(1,:)=[];
        SM_mat(1,:)=[];
        T_str=2.^(0:10);
        y_name='Weight magnitude';
    end

    subplot(2,2,kk);
    imagesc(correlation_mat); 
    if kk==4
        hc=colorbar('SouthOutside');
        cpos = get(hc,'position');
        cpos=cpos +[-0.25 -0.2 0 0];
        set(hc,'position',cpos);
    end
        set(gca, 'CLim', [0 1]); 
    if or(kk==1,kk==2)
        set(gca,'XTick', []);
    else
        set(gca,'XTick', 1:L_T, 'XTickLabel',T_str)
        xlabel('T [100 sec]','fontsize',labelfontsize); 
    end
  
    set(gca,'YTick', 1:length(y_str), 'YTickLabel',y_str);
    hy=ylabel(y_name);
    if kk==3         
         ypos = get(hy,'position');
        ypos=ypos +[-1.1 0 0];
        set(hy,'position',ypos)
    elseif kk==1
         ypos = get(hy,'position');
        ypos=ypos +[0.2 0 0];
        set(hy,'position',ypos)
    elseif kk==2
         ypos = get(hy,'position');
        ypos=ypos +[-0.3 0 0];
        set(hy,'position',ypos);
    end

    letter=['(' char(kk-1+'A') ')'] ;
    title(letter,'color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontweight','bold')

%     title('C','fontsize',labelfontsize);
end
% 
target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript\FinalProduction';
Export2Folder(['Fig7.tif'],target_folder) 