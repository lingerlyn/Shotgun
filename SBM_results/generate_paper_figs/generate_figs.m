    
set(0,'DefaultTextInterpreter', 'latex');
    
K=5; %width of subplots
x_ticks={'R','C','Z','S'};
fontsize=12;
fontsize2=1.5*fontsize;

L=2;
subplot = @(m,n,p) subtightplot (m, n, p, [0.07 0.03], .05,.04);
%% DIST DEP

load('regular_sbm_results','allEWs','W');
regEW=allEWs{1};
ddEW=allEWs{end};

   figure;
    ii=1;

    mi=min(W(:));ma=max(W(:));
    subplot(L+1,K,[1 2])
    imagesc(W,[mi ma]); h=colorbar;

    ylabel('True W','fontsize',fontsize2)
    set(h, 'ylim', [mi ma])
    subplot(L+1,K,[3 4])
    scatter(W(:),W(:),'b.'); box on;
    xlabel('W','fontsize',fontsize2)
    ylabel('W','fontsize',fontsize2)
    axis([mi ma mi ma])
    subplot(L+1,K,5)
    y_ticks=[1,1,1,1];   
    bar(y_ticks);    
    set(gca, 'XTickLabel', x_ticks,'fontsize',fontsize);
    
    %first estimate
    subplot(L+1,K,[K*ii+3 K*ii+4])
    A_ind=linspace(mi,ma,100);
    plot(A_ind,A_ind,'r-');
    hold all
    scatter(W(:),regEW(:),'b.')
    axis([mi ma mi ma])
    hold off
    xlabel('W','fontsize',fontsize2)
    ylabel('$\hat{W}$','fontsize',fontsize2)

    
   subplot(L+1,K,K*ii+5)
     [R,correlation, zero_matching,sign_matching] = GetWeightsErrors( W,regEW );

    bar( [R,correlation, zero_matching,sign_matching] );    
    ylim([0 1])
    set(gca, 'XTickLabel', x_ticks,'fontsize',fontsize);

    
    subplot(L+1,K,K*ii+[1 2])    
    imagesc(regEW,[mi ma]); h=colorbar;
    set(h, 'ylim', [mi ma])
    ylabel('Lasso','fontsize',fontsize2)
    
    %%% second estimate ii=2
    ii=2;
    subplot(L+1,K,K*ii+[3 4])
    A_ind=linspace(mi,ma,100);
    plot(A_ind,A_ind,'r-');
    hold all
    scatter(W(:),ddEW(:),'b.')
    axis([mi ma mi ma])
    hold off
    xlabel('W','fontsize',fontsize2)
    ylabel('$\hat{W}$','fontsize',fontsize2)

    
   subplot(L+1,K,K*ii+5)
     [R,correlation, zero_matching,sign_matching] = GetWeightsErrors( W,ddEW );

    bar( [R,correlation, zero_matching,sign_matching] );    
    ylim([0 1])
    set(gca, 'XTickLabel', x_ticks,'fontsize',fontsize);

    
    subplot(L+1,K,K*ii+[1 2])    
    imagesc(ddEW,[mi ma]); h=colorbar;
    set(h, 'ylim', [mi ma])
    ylabel('Lasso w. dist-dep','fontsize',fontsize2)
    %% REGULAR
    
load('regular_sbm_results','allEWs','W');
regEW=allEWs{1};
sbmEW=allEWs{end};

   figure;
    ii=1;

    mi=min(W(:));ma=max(W(:));
    subplot(L+1,K,[1 2])
    imagesc(W,[mi ma]); h=colorbar;

    ylabel('True W','fontsize',fontsize2)
    set(h, 'ylim', [mi ma])
    subplot(L+1,K,[3 4])
    scatter(W(:),W(:),'b.'); box on;
    xlabel('W','fontsize',fontsize2)
    ylabel('W','fontsize',fontsize2)
    axis([mi ma mi ma])
    subplot(L+1,K,5)
    y_ticks=[1,1,1,1];   
    bar(y_ticks);    
    set(gca, 'XTickLabel', x_ticks,'fontsize',fontsize);
    
    %first estimate
    subplot(L+1,K,[K*ii+3 K*ii+4])
    A_ind=linspace(mi,ma,100);
    plot(A_ind,A_ind,'r-');
    hold all
    scatter(W(:),regEW(:),'b.')
    axis([mi ma mi ma])
    hold off
    xlabel('W','fontsize',fontsize2)
    ylabel('$\hat{W}$','fontsize',fontsize2)

    
   subplot(L+1,K,K*ii+5)
     [R,correlation, zero_matching,sign_matching] = GetWeightsErrors( W,regEW );

    bar( [R,correlation, zero_matching,sign_matching] );    
    ylim([0 1])
    set(gca, 'XTickLabel', x_ticks,'fontsize',fontsize);

    
    subplot(L+1,K,K*ii+[1 2])    
    imagesc(regEW,[mi ma]); h=colorbar;
    set(h, 'ylim', [mi ma])
    ylabel('Lasso','fontsize',fontsize2)
    
    %%% second estimate ii=2
    ii=2;
    subplot(L+1,K,K*ii+[3 4])
    A_ind=linspace(mi,ma,100);
    plot(A_ind,A_ind,'r-');
    hold all
    scatter(W(:),sbmEW(:),'b.')
    axis([mi ma mi ma])
    hold off
    xlabel('W','fontsize',fontsize2)
    ylabel('$\hat{W}$','fontsize',fontsize2)

    
   subplot(L+1,K,K*ii+5)
     [R,correlation, zero_matching,sign_matching] = GetWeightsErrors( W,sbmEW );

    bar( [R,correlation, zero_matching,sign_matching] );    
    ylim([0 1])
    set(gca, 'XTickLabel', x_ticks,'fontsize',fontsize);

    
    subplot(L+1,K,K*ii+[1 2])    
    imagesc(sbmEW,[mi ma]); h=colorbar;
    set(h, 'ylim', [mi ma])
    ylabel('Lasso w. dist-dep','fontsize',fontsize2)