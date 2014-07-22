mi=min(W(:));ma=max(W(:));
%% Plot
figure
subplot(2,2,1); imagesc(Cxx); colorbar;
subplot(2,2,2); imagesc(Cxy); colorbar;
subplot(2,2,3); imagesc(EW,[mi ma]); h=colorbar;
set(h, 'ylim', [mi ma])
subplot(2,2,4); imagesc(W,[mi ma]); h=colorbar;
set(h, 'ylim', [mi ma])


% figure
% subplot(2,1,1); imagesc(sign(EW)); colorbar;
% subplot(2,1,2); imagesc(sign(W)); colorbar;

%%
A_ind=linspace(mi,ma,100);
figure
plot(A_ind,A_ind);
hold all
scatter(W(:),EW(:))
hold all
scatter(W(:),EW2(:))
legend('x=y','Direct','VB')
xlabel('True weights')
ylabel('Estimated weights')
[MSE,correlation,SIGN_ERROR] = GetWeightsErrors( W,EW )
[MSE2,correlation2,SIGN_ERROR2] = GetWeightsErrors( W,EW2 )

title({[' Direct corr =' num2str(correlation) ', VB corr =' num2str(correlation2)]; ...
     [' Direct MSE =' num2str(MSE) ', VB MSE =' num2str(MSE2)]; ...
     [' Direct SE =' num2str(SIGN_ERROR) ', VB SE =' num2str(SIGN_ERROR2) ]});
hold all
%%
% [~, ind]=max(EW(:));
% A_ind=linspace(-W(ind),W(ind),100); 
% A0=A_ind(end);
% y=A_ind/A0;
% f_A=log((1+y)./(1-y));
% hold all
% plot(A_ind,f_A);



