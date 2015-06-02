close all
clear all
clc

addpath('EstimateConnectivity')
addpath('Results/CommonInputFigure')
addpath('Plotting')
SetDefaultGraphicSettings(0)

K=16;
title_pos=[-0.02 0.95];
title_font=12;
set(0,'DefaultAxesFontSize',11)
h=figure(1)
% set(h,'units','normalized','outerposition',[  0.5130    0.1519    0.4042    0.6102])
units = 'centimeters';
set(gcf, 'PaperUnits', units,'Units', units)           
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperPosition',[4.5073   10.3273   13    9])
set(gcf, 'Position',[4.5073   10.3273   13    9])


subplot = @(m,n,p) subtightplot (m, n, p, [0.07 0.05], [0.12 0.04], [0.08 0.01]);

blue_red_map= darkb2r(-1,1);
colormap(blue_red_map)

load('Run_N=50_obs=1_T=5000000')
N=params.connectivity.N;

b=(V*Cxx)^(-1);
a=Cxy;
rates=V(eye(N)>0.5);
EW_obs=a'*b;
[amp, ~]=logistic_ELL(rates,EW_obs,Cxx,Cxy);
% amp=1
EW_obs=diag(amp)*EW_obs;
mi1=min(W(:));ma1=max(W(:));
if mi1>-ma1
    mi1=-ma1;
elseif ma1<-mi1
    ma1=-mi1;
end

subplot(2,2,1)
imagesc(W,[mi1 ma1])
colorbar
% title('true weights')
title('(A)', 'Units', 'normalized', 'Position', title_pos, 'HorizontalAlignment', 'right','fontweight','bold','fontsize',title_font) 
set(gca,'DataAspectRatio',[1 1 1])
 ylabel('neuron #');

x=0.5;y=0.5;w=16;h=16;
rectangle('Position',[x,y,w,h],'LineWidth',2,'edgecolor','w')

subplot(2,2,2)
temp=W(1:K,1:K);
mi=min(temp(:));ma=max(temp(:));
if mi>-ma
    mi=-ma;
elseif ma<-mi
    ma=-mi;
end
imagesc(temp,[mi ma])
colorbar
% title('true weights - subset')
title('(B)', 'Units', 'normalized', 'Position',title_pos, 'HorizontalAlignment', 'right','fontweight','bold','fontsize',title_font) 
set(gca,'DataAspectRatio',[1 1 1])



b=(V(1:K,1:K)*Cxx(1:K,1:K))^(-1);
a=Cxy(1:K,1:K);
rates=V(eye(N)>0.5);
EW_unobs=a'*b;
[amp, Ebias2]=logistic_ELL(rates(1:K),EW_unobs,Cxx(1:K,1:K),Cxy(1:K,1:K));
% amp=1;
EW_unobs=diag(amp)*EW_unobs;



subplot(2,2,4)
imagesc(EW_unobs,[mi ma])
% title('infered weights - subset, with hidden neurons')
title('(D)', 'Units', 'normalized', 'Position', title_pos, 'HorizontalAlignment', 'right','fontweight','bold','fontsize',title_font)
set(gca,'DataAspectRatio',[1 1 1])
 xlabel('neuron #');
colorbar

load('Run_N=50_obs=0.3_T=5000000.mat')
b=(V*Cxx)^(-1);
a=Cxy;
rates=V(eye(N)>0.5);
EW_obs=a'*b;
[amp, Ebias2]=logistic_ELL(rates,EW_obs,Cxx,Cxy);
% amp=1
EW_obs=diag(amp)*EW_obs;

subplot(2,2,3)
imagesc(EW_obs(1:K,1:K),[mi ma])
% title('infered weights - subset, with shotgun')
title('(C)', 'Units', 'normalized', 'Position',title_pos, 'HorizontalAlignment', 'right','fontweight','bold','fontsize',title_font) 
set(gca,'DataAspectRatio',[1 1 1])
 ylabel('neuron #');
  xlabel('neuron #');
colorbar

%%
target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript\Revision2'
Export2Folder(['Fig3.eps'],target_folder) 

%%
mean(EW_obs(~W(1:K,1:K)))
mean(EW_unobs(~W(1:K,1:K)))
%%
std(EW_obs(~W(1:K,1:K)))
std(EW_unobs(~W(1:K,1:K)))
%%
% t=1:50;
% subplot(2,1,1)
% plot(t,U(1,t),t,U(2,t))
% subplot(2,1,2)
% plot(t,spikes(1,t),'xb',t,spikes(2,t),'og')