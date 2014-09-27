K=2;
% U=W*spikes;
% corr(U(1:K,:)',U(1:K,:)')
% corr(spikes(1:K,:)',spikes(1:K,:)')

% corr(U(1:K,1:end-1)',U(1:K,2:end)')
% corr(spikes(1:K,1:end-1)',spikes(1:K,2:end)')

Cxx(1:K,1:K);
b=(V(1:K,1:K)*Cxx(1:K,1:K))^(-1);
a=Cxy(1:K,1:K);
temp=W(1:K,1:K);
EW3_obs=EW3(1:K,1:K);
figure
mi=min(temp(:));ma=max(temp(:));
subplot(3,1,1)
imagesc(W(1:K,1:K),[mi ma])
colorbar
title('weight matrix')
subplot(3,1,2)
imagesc(EW3_obs,[mi ma])
title('common input observed')
colorbar
EW3_unobs=a'*b;
subplot(3,1,3)
imagesc(EW3_unobs,[mi ma])
title('common input unobserved')
colorbar
%%
mean(EW3_unobs(~W(1:K,1:K)))
mean(EW3_obs(~W(1:K,1:K)))
%%
std(EW3_unobs(~W(1:K,1:K)))
std(EW3_obs(~W(1:K,1:K)))
%%
% t=1:50;
% subplot(2,1,1)
% plot(t,U(1,t),t,U(2,t))
% subplot(2,1,2)
% plot(t,spikes(1,t),'xb',t,spikes(2,t),'og')