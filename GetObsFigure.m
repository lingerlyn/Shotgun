clear all
% close all
clc

% Generate all figures for the paper
addpath('GenerateSpikes');
addpath('Plotting')
addpath('Misc')
SetDefaultGraphicSettings(1)

subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.03], [0.2 0.1], [0.1 0.1]);
% subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.08], [0.07 0.04], [0.05 0.01]);

title_pos=[-0.11,0.85];
title_pos2=[-0.25,0.86];
figure(1)
% set(h,'units','normalized','outerposition',[0 0 0.6 1])
units = 'centimeters';
set(gcf, 'PaperUnits', units,'Units', units)           
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperPosition',[-0.0219   -0.2626   15   22])
set(gcf, 'Position',[-0.0219   -0.2626   15   22])


N=10;
T=1e5*N;
t_show=1:3e3;
dt=0.01; %100 Hz frame rate
tt=t_show*dt;
sample_ratio=0.2;
N_stim=0;
sample_type_set={'continuous','fixed_subset','spatially_random','prob','double_continuous','fully_random'};
obs_duration=100;
t_start=1;

seed_sample=1;
L=100; %number of histogram bins
M=5; %show M cases
a=M; b=5; %subplot grid

for kk=1:M
switch kk
    case 1
        sample_type=sample_type_set{2};
        name='Fixed';
    case 2
        sample_type=sample_type_set{1};
        name='Serial';
    case 3
        sample_type=sample_type_set{6};
        name='Fully random';        
    case 4
        sample_type=sample_type_set{3};
        name='Random blocks';
    case 5
        sample_type=sample_type_set{5};
        name='Double serial';
end 

observations=SampleSpikes(N,T,sample_ratio,sample_type,obs_duration,N_stim,seed_sample,t_start);
% inputs:
subplot(a,b,(kk-1)*b+[1 2 3])
imagesc(tt,[],observations(:,t_show))
% title('(A)', 'Units', 'normalized', 'Position', [title_horz title_height], 'HorizontalAlignment', 'right','fontweight','bold','fontsize',title_font) 
if kk==M
%     xlabel('time [sec]','fontweight','bold');
        xlabel('time [sec]')
else
    set(gca,'xtick',[]);
end
% ylabel(name);
ylabel('neuron #');

colormap('gray')
freezeColors
letter=['(' char( 2*(kk-1)+'A') ')'] ;
title(letter,'color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos,'fontweight','bold')

subplot(a,b,(kk-1)*b+[4 5])
sample_type=sample_type_set{3};
% observations=SampleSpikes(N,T,sample_ratio,sample_type,obs_duration,N_stim,seed_sample);
XYn=observations(:,1:(end-1))*(observations(:,2:end))'/T;
ma=max(XYn(:));
imagesc(XYn,[0 ma]);
set(gca,'DataAspectRatio',[1 1 1])
ylabel('neuron #');
caxis([0 0.05]);
% ['random scanning' -   ; 'all spike pairs observed (good)']
letter=['(' char( 2*(kk-1)+1+'A') ')'] ;
title(letter,'color', 'k', 'Units', 'normalized', 'interpreter','none','position',title_pos2,'fontweight','bold')

colormap(paruly)

freezeColors
if kk==1
    hc=colorbar('location','northoutside');
    cpos = get(hc,'position');
    cpos=cpos +[-0.045 0.04 0.085 0];
    set(hc,'position',cpos);
%     allh=get(gcf, 'Children');
%     hc=allh(ismember(get(allh,'Tag'),'Colorbar'));
%         hc=cbfreeze(hc)
end

if kk==M
    set(hc,'XTick',[0, 0.025, 0.05]);    
%     colorbar_labels = get(hc,'XTickLabel');
    colorbar_labels=['  0  '; '0.025'; ' 0.05'];
    colorbar_labels=[colorbar_labels, [repmat(' ',1,size(colorbar_labels,1)-1), '<']'];
    set(hc,'XTickLabel',colorbar_labels);
end


if kk==M
    xlabel('neuron #');
else
    set(gca,'xtick',[]);
end
end

%%
target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript\Revision2';
figure_name='Fig1.eps';
    
set(gcf, 'Color', 'w'); 
set(findall(gcf,'type','axe'),'box','on')
% set(findall(gcf,'type','text'),'fontWeight','bold')
%     set(gca,'fontWeight','bold')
print(gcf,figure_name, '-painters','-depsc2')

movefile(figure_name, fullfile(target_folder,figure_name));
