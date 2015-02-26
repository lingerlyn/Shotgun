close all

dt=0.01
T=size(Y,2);
N=size(Y,1);
tt=(1:T)*dt;
ii=6;
figure(100)
plot(tt,Y(ii,:)/max(Y(ii,:)),'k')
xlim([0 20])
hold all
plot(tt,spikes(ii,:),'ob',tt,true_spikes(ii,:),'.r')
ylim([0 1.1]);
%%
figure(200)
subplot(3,1,1)
imagesc(true_spikes)
subplot(3,1,2)
imagesc(Y)
subplot(3,1,3)
imagesc(spikes)

%%
figure(300)
correlation=zeros(N,1);
for ii=1:N
correlation(ii)=corr(spikes(ii,:)',true_spikes(ii,:)');
end
plot(correlation)
find(isnan(correlation))
