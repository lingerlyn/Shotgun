% close all

dt=0.01;
T=size(Y,2);
N=size(Y,1);
tt=(1:T)*dt;
ii=2;
figure(100)
plot(tt,Y(ii,:)/max(Y(ii,:)),'k')
% xlim([0 tt(end)])
xlim([0 20])

hold all
figure(20)
plot(tt,spikes(ii,:),'ob',tt,true_spikes(ii,:),'.r')
ylim([0 1.1]);

%% spike correlation 
figure(300)
Spike_rec_correlation=zeros(N,1);
for ii=1:N
    Spike_rec_correlation(ii)=corr(spikes(ii,:)',true_spikes(ii,:)');
end
plot(Spike_rec_correlation)
find(isnan(Spike_rec_correlation))

%% binned  spike correlation 
bin_min=1;bin_max=30;
DF = bin_min:bin_max;
CR = zeros(N,length(DF));
KL = CR;
for nn = 1:N
    [CR(nn,:),KL(nn,:)] = ca_metrics(squeeze(spikes(nn,:)'),true_spikes(nn,:)',DF);
end

timebins=DF*dt;
figure(999)
subplot(2,1,1)
plot(timebins,CR)
ylabel('Correlation')
xlabel('Time bin size [sec]')   
xlim([bin_min bin_max]*dt)
subplot(2,1,2)
plot(timebins,KL)
ylabel('-KL divergence')
xlabel('Time bin size [sec]')    
xlim([bin_min bin_max]*dt)

%%
% figure(200)
% subplot(3,1,1)
% imagesc(true_spikes)
% subplot(3,1,2)
% imagesc(spikes)
% subplot(3,1,3)
% imagesc(Y)

%% 
% figure(201)
% hold off
% for kk=1:N
%     subplot(2,1,1)
%     plot(relative_std_cell{kk})
%     hold all
%     subplot(2,1,2)
%     plot(diff(relative_std_cell{kk}))
%     hold all
% end


%% Compare covariances qualites

mY=full(sum(true_spikes,2));
XX=true_spikes*true_spikes';
XY=true_spikes(:,1:(end-1))*(true_spikes(:,2:end))';

rates_true=full(mY./(mYn+eps)); %estimate the mean firing rates
        
Cxx_true=full(XX./(XXn+eps))-rates_true*rates_true'; %estimate the covariance (not including stim terms for now)
Cxy_true=full(XY./(XYn+eps))-rates_true*rates_true'; %estimate cross-covariance

%% Attempt to correct covariances using missing/added spikes error model
d_noise=full(mean(spikes.*true_spikes,2)./mean(true_spikes,2)); % missing spikes - later should be estimated from data
e_noise=full(mean(spikes.*(1-true_spikes),2)./mean(1-true_spikes,2)); 
mY=full(sum(sampled_spikes_obs,2));
rates_adjusted=(mY./(mYn+eps)-e_noise)./(d_noise-e_noise); %estimate the mean firing rates
C_mat=1-bsxfun(@plus,rates_adjusted',rates_adjusted);

XX=sampled_spikes_obs*sampled_spikes_obs';
XY=sampled_spikes_obs(:,1:(end-1))*(sampled_spikes_obs(:,2:end))';
E_mat=e_noise*e_noise';
D_mat=d_noise*d_noise';
Cxy_adjusted=(full(XY./(XYn+eps))-C_mat.*E_mat)./(E_mat+D_mat)-rates_adjusted*rates_adjusted'; %estimate cross-covariance

E_mat(eye(N)>0.5)=e_noise;
d_mat(eye(N)>0.5)=d_noise;
Cxx_adjusted=(full(XX./(XXn+eps))-C_mat.*E_mat)./(E_mat+D_mat)-rates_adjusted*rates_adjusted'; 


%% First Attempt to correct covariances
% c=1./Spike_rec_correlation; %correction factor
% mY=full(sum(sampled_spikes_obs,2)).*c;
% XX=sampled_spikes_obs*sampled_spikes_obs';
% temp=XX(eye(N_obs)>0.5).*c;
% XX=XX.*(c*c');
% XX(eye(N_obs)>0.5)=temp;
% XY=sampled_spikes_obs(:,1:(end-1))*(sampled_spikes_obs(:,2:end))';
% XY=XY.*(c*c');
% rates_adjusted=mY./(mYn+eps); %estimate the mean firing rates
% Cxx_adjusted=full(XX./(XXn+eps))-rates_adjusted*rates_adjusted'; 
% Cxy_adjusted=full(XY./(XYn+eps))-rates_adjusted*rates_adjusted'; %estimate cross-covariance

%% plot
figure(522)
a=2;b=4;
subplot(a,b,1)
mi_xx=min(Cxx(:));ma_xx=max(Cxx(:));
imagesc(Cxx,[mi_xx,ma_xx])

subplot(a,b,2)
imagesc(Cxx_true,[mi_xx,ma_xx])

subplot(a,b,3)
imagesc(Cxx_adjusted,[mi_xx,ma_xx])

subplot(a,b,4)
mi=min([Cxx_true(:)])*1.2;
ma=max([Cxx_true(:)])*1.2;
A_ind=linspace(mi,ma,100);
plot(A_ind,A_ind,'k--','linewidth',1);
hold all
plot(Cxx_true(:),Cxx(:),'.',Cxx_true(:),Cxx_adjusted(:),'.')

subplot(a,b,5)
mi_xy=min(Cxy(:));ma_xy=max(Cxy(:));
imagesc(Cxy,[mi_xy,ma_xy])

subplot(a,b,6)
imagesc(Cxy_true,[mi_xy,ma_xy])

subplot(a,b,7)
imagesc(Cxy_adjusted,[mi_xy,ma_xy])

subplot(a,b,8)
mi=min([Cxy_true(:)])*1.2;
ma=max([Cxy_true(:)])*1.2;
A_ind=linspace(mi,ma,100);
plot(A_ind,A_ind,'k--','linewidth',1);
hold all
plot(Cxy_true(:),Cxy(:),'.',Cxy_true(:),Cxy_adjusted(:),'.')