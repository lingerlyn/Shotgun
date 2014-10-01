%% known mu
for oo=1:nobs
figure(100); hold on;
h1=plot(squeeze(omp_corrs(1,oo,:)),'b'); plot(squeeze(lasso_ell_corr(1,oo)),'ob');
h2=plot(squeeze(omp_corrs(2,oo,:)),'g'); plot(squeeze(lasso_ell_corr(2,oo)),'og');
h3=plot(squeeze(omp_corrs(3,oo,:)),'c'); plot(squeeze(lasso_ell_corr(3,oo)),'oc');
% plot(squeeze(omp_ell_corrs(1,oo,:)),'--b');
% plot(squeeze(omp_ell_corrs(2,oo,:)),'--g');
% plot(squeeze(omp_ell_corrs(3,oo,:)),'--c');

title(['obs frac = ' num2str(obs(oo))]); xlabel('lambda index'); ylabel('correlation')
legend([h1,h2,h3],'T1e4','T5e4','T1e5')
xlim([.75,nls]);

saveas(100,['oo' num2str(oo) '_corr'],'fig');
saveas(100,['oo' num2str(oo) '_corr'],'png');
close;

figure(101); hold on;
h1=plot(squeeze(omp_mses(1,oo,:)),'b'); plot(squeeze(lasso_ell_mse(1,oo)),'ob');
h2=plot(squeeze(omp_mses(2,oo,:)),'g'); plot(squeeze(lasso_ell_mse(2,oo)),'og');
h3=plot(squeeze(omp_mses(3,oo,:)),'c'); plot(squeeze(lasso_ell_mse(3,oo)),'oc');
% plot(squeeze(omp_ell_mses(1,oo,:)),'--b');
% plot(squeeze(omp_ell_mses(2,oo,:)),'--g');
% plot(squeeze(omp_ell_mses(3,oo,:)),'--c');

title(['obs frac = ' num2str(obs(oo))]); xlabel('lambda index'); ylabel('relative mse');
legend([h1,h2,h3],'T1e4','T5e4','T1e5')
xlim([.75,nls]);

saveas(101,['oo' num2str(oo) '_mse'],'fig');
saveas(101,['oo' num2str(oo) '_mse'],'png');
close;


%%%%% daniel's figure %%%%%%%%%
for tt=1:nts

ll=3;

EW=omp_ests(:,:,tt,oo,ll);
EW2=lasso_ell_ests(:,:,tt,oo);

figure(103);
subplot(2,3,1);

imagesc(EW); colorbar; caxis([min(W(:)),max(W(:))]);

subplot(2,3,2); 
% ymin=min([EW(:);EW2(:)]);
% ymax=max([EW(:);EW2(:)]);

plot(W(:),EW(:),'o');
% ylim([ymin,ymax]);
% xlim([min(W(:)),max(W(:))]);
xlabel('True W'); ylabel('Estimated W');
title(['omp; l = ' num2str(ls(ll)) ])

subplot(2,3,3); [R,C,Z,S] = GetWeightsErrors( W,EW );
bar(1:4,[R,C,Z,S]);
set(gca,'XTickLabel',{'R', 'C', 'Z', 'S'})

subplot(2,3,4);
imagesc(EW2); colorbar; caxis([min(W(:)),max(W(:))]);

subplot(2,3,5);
plot(W(:),EW2(:),'o');
% ylim([ymin,ymax]);
% xlim([min(W(:)),max(W(:))]);
xlabel('True W'); ylabel('Estimated W');
title('lasso + ELL')

subplot(2,3,6); [R,C,Z,S] = GetWeightsErrors( W,EW2 );
bar(1:4,[R,C,Z,S]);
set(gca,'XTickLabel',{'R', 'C', 'Z', 'S'})

suptitle(['T = ' num2str(ts(tt)) ';obs ' num2str(obs(oo))]);
saveas(103,['tt' num2str(tt) 'oo' num2str(oo) '_comp'],'fig');
saveas(103,['tt' num2str(tt) 'oo' num2str(oo) '_comp'],'png');
end

end

%% unknown mu looping

oo=4;
tt=2;
ll=3;

figure; hold on;
plot(0:niter-1,loop_mses)
plot(0,omp_mses(tt,oo,ll),'*r')
plot(0,lasso_ell_mse(tt,oo),'o');
title('mse');
xlabel('iteration');
legend 'using inferred M' 'using real M' 'lasso solution'

figure; hold on;
plot(0:niter-1,loop_corrs)
plot(0,omp_corrs(tt,oo,ll),'*r')
plot(0,lasso_ell_corr(tt,oo),'o');
title('correlation')
xlabel('iteration')
legend 'using inferred M' 'using real M' 'lasso solution'

figure;
plot(M_mses);
xlabel('iteration')
title('Mean matrix mses');

figure; plot(M_corrs);
xlabel('iteration')
title('Mean matrix correlations')


%%% daniel figure
%%
EW=allEWs{end};
EW2=lasso_ell_ests(:,:,tt,oo);

figure(104);
subplot(2,3,1);

imagesc(EW); colorbar; caxis([min(W(:)),max(W(:))]);

subplot(2,3,2); 
% ymin=min([EW(:);EW2(:)]);
% ymax=max([EW(:);EW2(:)]);

plot(W(:),EW(:),'o');
% ylim([ymin,ymax]);
% xlim([min(W(:)),max(W(:))]);
xlabel('True W'); ylabel('Estimated W');
title(['omp; l = ' num2str(ls(ll)) ])

subplot(2,3,3); [R,C,Z,S] = GetWeightsErrors( W,EW );
bar(1:4,[R,C,Z,S]);
set(gca,'XTickLabel',{'R', 'C', 'Z', 'S'})

subplot(2,3,4);
imagesc(EW2); colorbar; caxis([min(W(:)),max(W(:))]);

subplot(2,3,5);
plot(W(:),EW2(:),'o');
% ylim([ymin,ymax]);
% xlim([min(W(:)),max(W(:))]);
xlabel('True W'); ylabel('Estimated W');
title('lasso + ELL')

subplot(2,3,6); [R,C,Z,S] = GetWeightsErrors( W,EW2 );
bar(1:4,[R,C,Z,S]);
set(gca,'XTickLabel',{'R', 'C', 'Z', 'S'})

suptitle(['T = ' num2str(ts(tt)) ';obs ' num2str(obs(oo))]);