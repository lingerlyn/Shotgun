clear all
close all
clc

SetDefaultGraphicSettings(1)

h=figure(1)
% set(h,'Units','centimeters','Position',[2 5 20 10])
units = 'centimeters';
set(gcf, 'PaperUnits', units,'Units', units)           
set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'PaperPosition',[-0.0219   -0.2626   13   8])
set(gcf, 'Position',[-0.0219   -0.2626   13   8])

x=-5:0.001:5;
c=sqrt(pi/8);


subplot = @(a,b,p) subtightplot (a,b, p, [0.06 0.08], [0.2 0.06], [0.08 0.05]);
% Plot Log approximation
subplot(1,2,1)
s_array=[0.25,0.5,1,2,4,8];
lgnd_str={};

for s=s_array
mu=-5:0.01:5;

W=bsxfun(@times,s,randn(50000,1));
X=bsxfun(@plus,mu,W);
y=mean(log((1+exp(X))),1);
z=sqrt(1+(pi*s.^2)/8).*log((1+exp(mu./(sqrt(1+(pi*s.^2)/8)))));
    colorOrder = get(gca, 'ColorOrder');
    current_color=colorOrder(mod(length(get(gca, 'Children')), size(colorOrder, 1))+1, :);
% plot(mu,y,'-',mu,z,'--','color',current_color)
plot(mu,y./z,'-','color',current_color)
ylim([0, 1.1])
xlabel('\mu')
hold all
% lgnd_str{end+1}=['$\sigma$=' num2str(s)]; %#ok
end
% legend(lgnd_str,'location','southeast');
set(gca,'xlim',[mu(1),mu(end)],'xtick',mu(1):mu(end))
hold off


% Plot Sigmoid approximation
subplot(1,2,2)
for s=s_array
mu=-5:0.01:5;

W=bsxfun(@times,s,randn(50000,1));
X=bsxfun(@plus,mu,W);
y=mean(1./((1+exp(-X))),1);
z=1./(1+exp(-mu./sqrt(1+(pi/8)*s.^2)));
    colorOrder = get(gca, 'ColorOrder');
    current_color=colorOrder(mod(length(get(gca, 'Children')), size(colorOrder, 1))+1, :);
% plot(mu,y,'-',mu,z,'--','color',current_color)
plot(mu,y./z,'-','color',current_color)
xlabel('\mu')
hold all
lgnd_str{end+1}=['\sigma=' num2str(s)]; %#ok
end
legend(lgnd_str,'location','southeast');
set(gca,'xlim',[mu(1),mu(end)],'xtick',mu(1):mu(end))
ylim([0, 1.1])
hold off

target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript\Revision2';
Export2Folder(['WangApprox.eps'],target_folder) 

%%
% for s=s_array
% mu=-5:0.01:5;
% 
% W=bsxfun(@times,s,randn(50000,1));
% X=bsxfun(@plus,mu,W);
% y=mean(1./((1+exp(-X))),1);
% z=sigmoid_int(mu,s);
%     colorOrder = get(gca, 'ColorOrder');
%     current_color=colorOrder(mod(length(get(gca, 'Children')), size(colorOrder, 1))+1, :);
% % plot(mu,y,'-',mu,z,'--','color',current_color)
% plot(mu,y./z,'-','color',current_color)
% xlabel('\mu')
% hold all
% lgnd_str{end+1}=['\sigma=' num2str(s)]; %#ok
% end
% legend(lgnd_str,'location','southeast');
% set(gca,'xlim',[mu(1),mu(end)],'xtick',mu(1):mu(end))
% ylim([0, 1.1])
% hold off



% %% Other OLD stuff
% 
% 
% % c=sqrt(2*pi)*exp(b^2/2)/(2+exp(b)+exp(-b));
% % y=normcdf(c*x,0,1);
% % c=0.36;
% % y=(x>0)-((x>0)-normcdf(c*x,0,1)).*(exp(-abs(c*x)));
% % y=(x>0).*(1-exp(-x))+(x<0).*exp(x);
% b=2;
% % y=(x>b).*(1-exp(-x))+normcdf(c*x,0,1).*(x>=-b).*(x<=b)+(x<-b).*exp(x);
% y=(x>b).*(1-exp(-x))+(0.5+0.25*x-x.^3/48+x.^5/480).*(x>=-b).*(x<=b)+(x<-b).*exp(x);
% z=1./(1+exp(-x));
% 
% % subplot(1,2,1)
% plot(x,y,'-b',x,z,':r')
% title('linear scale')
% %%
% x=-10:0.01:2;
% y=(x>0)-((x>0)-normcdf(c*x,0,1)).*(exp(-abs(c*x)));
% z=1./(1+exp(-x));
% subplot(1,2,2)
% semilogy(x,y,'-b',x,z,':r')
% xlim([min(x) max(x)])
% % title('log scale')
%  
% % target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript';
% % Export2Folder(['WangApprox.eps'],target_folder) 
% 
% %% version 2
% clear all
% close all
% clc
% 
% set(0,'DefaultTextInterpreter', 'latex');
% set(0,'DefaultAxesFontSize',17)
% set(0,'defaultlinelinewidth',2.5)
% 
% 
% h=figure
% set(h,'Units','centimeters','Position',[2 5 20 10])
% x=-5:0.001:5;
% c=sqrt(pi/8);
% % c=0.36;
% % c=sqrt(2*pi)*exp(b^2/2)/(2+exp(b)+exp(-b));
% y=normcdf(c*x,0,1);
% 
% % y=(x>0)-((x>0)-normcdf(c*x,0,1)).*(exp(-abs(c*x)));
% z=1./(1+exp(-x));
% 
% subplot(1,2,1)
% plot(x,y,'-b',x,z,':r')
% title('linear scale')
% subplot(1,2,2)
% semilogy(x,1-y,'-b',x,1-z,':r')
% title('log scale')
% 
% % target_folder='C:\Users\Daniel\Copy\Columbia\Research\Shotgun\Manuscript';
% % Export2Folder(['WangApprox.eps'],target_folder) 
% %%
% L=1e3;
% M=100;
% b=0.01;
% Sigma=1;
% y=zeros(M,1);
% z=y;
% mu=linspace(-10,10,M);
% for ii=1:M
%     x=mu(ii)+Sigma*randn(L,1);
%     z(ii)=mean(1./(1+exp(-x)));
%     y(ii)=sigmoid_int(mu(ii),Sigma);
% end
% plot(mu,z,mu,y);
% % semilogy(mu,z,mu,y);