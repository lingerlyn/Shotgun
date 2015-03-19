close all

dt=0.01
T=size(Y,2);
N=size(Y,1);
tt=(1:T)*dt;
ii=4;
figure(100)
plot(tt,Y(ii,:)/max(Y(ii,:)),'k')
xlim([0 20])
hold all
plot(tt,spikes(ii,:),'ob',tt,true_spikes(ii,:),'.r')
% ylim([0 1.1]);
%%
figure(200)
subplot(3,1,1)
imagesc(true_spikes)
subplot(3,1,2)
imagesc(spikes)
subplot(3,1,3)
imagesc(Y)

%% 
figure(201)
hold off
for kk=1:N
    subplot(2,1,1)
    plot(relative_std_cell{kk})
    hold all
    subplot(2,1,2)
    plot(diff(relative_std_cell{kk}))
    hold all
end

%%
figure(300)
Spike_rec_correlation=zeros(N,1);
for ii=1:N
    Spike_rec_correlation(ii)=corr(spikes(ii,:)',true_spikes(ii,:)');
end
plot(Spike_rec_correlation)
find(isnan(Spike_rec_correlation))
