[hist_real,bins1]=hist(diag(Cxx(ind_sample,ind_sample)),100);
[hist_sim,bins2]=hist(z(:),100);
plot(bins1,hist_real/sum(hist_real),bins2,hist_sim/sum(hist_sim))
ylabel('Prob Dist.'); xlabel('Cxx values')
legend('simulation','independent')

%% 
[hist_real,bins1]=hist(Cxx(ind_sample,ind_sample),100);
[hist_sim,bins2]=hist(z,100);
plot(bins1,hist_real/sum(hist_real),bins2,hist_sim/sum(hist_sim))
ylabel('Prob Dist.'); xlabel('Cxy values')
legend('simulation','independent')


%%
subplot(2,1,1)
hist(Cxy(W'==0),100)
subplot(2,1,2)
hist(Cxy(W'~=0),100)


%%
[hist_real,bins1]=hist(Cxy(W'==0),100);
[hist_sim,bins2]=hist(z,100);
plot(bins1,hist_real/sum(hist_real),bins2,hist_sim/sum(hist_sim))
ylabel('Prob Dist.'); xlabel('Cxy values')
legend('simulation','independent')

%%
[~,idx]=sort(W(:));
temp=Cxy;
tt=1:length(W(:));
plotyy(tt,temp(idx),tt,W(idx));
hold all

%%
[~,idx]=sort(Cxy_single(:));
plot(tt,Cxy_single(idx));
hold all

%%
plot(Cxy_single)
