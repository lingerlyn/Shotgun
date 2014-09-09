mi=min(W(:));ma=max(W(:));

%% Plot
figure(1)
subplot(2,2,1); imagesc(W,[mi ma]); h=colorbar;
title('True W')
set(h, 'ylim', [mi ma])
subplot(2,2,2); imagesc(EW,[mi ma]); h=colorbar;
title('EW1')
set(h, 'ylim', [mi ma])
subplot(2,2,3); imagesc(EW2,[mi ma]); h=colorbar;
title('EW2')
set(h, 'ylim', [mi ma])
subplot(2,2,4); imagesc(glassoEW,[mi ma]); h=colorbar;
title('glassoEW')
set(h, 'ylim', [mi ma])



% figure
% % subplot(2,1,1); imagesc(Cxx); colorbar;
% % subplot(2,1,2); imagesc(Cxy); colorbar;


% figure
% subplot(2,1,1); imagesc(sign(EW)); colorbar;
% subplot(2,1,2); imagesc(sign(W)); colorbar;

%%
A_ind=linspace(mi,ma,100);
figure(2)
plot(A_ind,A_ind);
hold all
scatter(W(:),EW(:),'.')
hold all
scatter(W(:),EW2(:),'.')
% hold all
% scatter(W(:),glassoEW(:),'.')
legend('x=y','EW','EW2','glassoEW')
xlabel('True weights')
ylabel('Estimated weights')
[R_squared,correlation,SE] = GetWeightsErrors( W,EW );
[R_squared2,correlation2,SE2] = GetWeightsErrors( W,EW2 );

title({[' EW corr =' num2str(correlation) ', EW2 corr =' num2str(correlation2)]; ...
     [' EW MSE =' num2str(R_squared) ', EW2 MSE =' num2str(R_squared2)]; ...
     [' EW SE =' num2str(SE) ', EW2 SE =' num2str(SE2) ]});
hold off

%%
mi=min(bias); ma=max(bias);
b_ind=linspace(mi,ma,100);
figure(3)
plot(b_ind,b_ind);
hold all
scatter(bias(:),Ebias(:))
hold all
scatter(bias(:),Ebias2(:))
legend('x=y','Ebias','Ebias2')
xlabel('True bias')
ylabel('Estimated bias')
[R_squared,correlation,SE] = GetWeightsErrors( bias,Ebias );
[R_squared2,correlation2,SE2] = GetWeightsErrors( bias,Ebias2 );

title({[' Eb corr =' num2str(correlation) ', Eb2 corr =' num2str(correlation2)]; ...
     [' Eb MSE =' num2str(R_squared) ', Eb2 MSE =' num2str(R_squared2)]; ...
     [' Eb SE =' num2str(SE) ', Eb2 SE =' num2str(SE2) ]});

hold off
%% Activity
figure(4)
subplot(2,3,[1 2])
imagesc(spikes);  colorbar
xlabel('time');
ylabel('Neurons');
firing_rate=mean(spikes(:));
title(['Spikes. Firing rate=' num2str(firing_rate)]);
subplot(3,3,3)
L=0:max(spikes(:));
[count,bins]=hist(spikes(:),L);  colorbar
stem(bins,count/(T*N));
title('S count')
subplot(2,3,[4 5])
U=bsxfun(@plus,W*spikes,bias);
imagesc(U); colorbar
xlabel('time');
ylabel('Neurons');
title(['U. E[exp(U)]=' num2str(mean(exp(U(:))))]);
subplot(3,3,6)
L=1000;
U_no_bias=W*spikes;
[count,bins]=hist(U_no_bias(:),L);  colorbar
% p0=exp(-spar*firing_rate*N); %probability for zero input
% if p0>1 
%     p0=1;
% end
mean_U=mean(mean(U_no_bias,2));
var_U=mean(var(U_no_bias,[],2));
pred_U_no_bias=exp(-(bins-mean_U).^2/(2*var_U));
pred_U_no_bias=pred_U_no_bias/sum(pred_U_no_bias);
plot(bins,count/(T*N),bins,pred_U_no_bias); %,bins,p0
title('U-b count')
subplot(3,3,9)
a=mean(exp(U),2);
b=exp(mean(U,2)+0.5*var(U,[],2));
c=exp(mean(U,2)).*(1+0.5*var(U,[],2));
plot(a,b,'.',a,c,'or',a,a,'-')
title('Expected loglikelihood Approximation check')
xlabel('E[exp(U)]')
ylabel('exp( \mu + 0.5\sigma^2 )')
hold off
%% Check Gaussianty
% figure
% var_U=var(U_no_bias,[],2);
% kar_U=mean((bsxfun(@plus,U_no_bias,-mean(U_no_bias,2))).^4,2);
% plot(kar_U,3*var_U.^2,'.',kar_U,kar_U,'-')