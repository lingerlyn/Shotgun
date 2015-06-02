clear all
% close all
clc


addpath('Misc')
addpath('EstimateConnectivity')
addpath('Results')

obs_array=[0.01];

for ii=1%1:length(obs_array)
    ii %#ok show iteration
    load(['Run_N=1100_obs=' num2str(obs_array(ii)) '_T=2000000_continuous_Cavity.mat'],'Cxx','Cxy','rates','W','bias','params','centers');
    % load('Run_N=1000_obs=1_T=2000000_Cavity','Cxx','Cxy','rates','W','params','centers');
    % load('Run_N=50_obs=1_T=100000_continuous_Cavity','Cxx','Cxy','rates','W','params','centers');

    [pen_diag,warm,est_type,~]=v2struct(params.conn_est_flags);
    pen_dist=0;
    N_stim=0; N=params.connectivity.N;
    est_spar=nnz(W(eye(N)<0.5))/(N^2-N);
    
    zeros_ind=find(W(end,:)==0);
    %% add inputs
    N_zero=round(N/est_spar); %additional zero inputs required
    ind_sample=zeros_ind(randi(length(zeros_ind),[N_zero,1])); %sample N neurons 

    %  Target input connectivity
    W_target=[zeros(1,N_zero), W(N,:)];

    % Input rates
    mock_rates=rates(ind_sample);
    rates_single=[mock_rates; rates];

    % Input correlations
    m=mock_rates-mock_rates.^2;
    s=sqrt(mock_rates.*(1-mock_rates)/params.spike_gen.T*params.spike_gen.sample_ratio);
    z=m+s.*randn(N_zero,1);
    Cxx_diag_single= [z; diag(Cxx)];

    % Input-output correlations
    m=Cxy(ind_sample,end)+mock_rates*rates(end);
    
    %indepedent model
    s=sqrt(m.*(1-m)/params.spike_gen.T*params.spike_gen.sample_ratio);
    z=Cxy(ind_sample,end)+s.*randn(N_zero,1);
    Cxy_single=[z; Cxy(:,end)];

    %%
    [EW,Ebias,quality,error_rates,lambda_path]=EstimateA_L1_logistic_cavity_SingleNeuron(Cxx_diag_single,Cxy_single,rates_single,est_spar,N_stim,pen_diag,pen_dist,warm,W_target,centers);     

    %% Plot
    addpath('Misc')
    addpath('Results')

    c=sum(abs(Cxx(:)))/sum(diag(Cxx));
    % c=1;
    mi=min(W(:));ma=max(W(:));
    A_ind=linspace(mi,ma,100);
    hold off
    plot(A_ind,A_ind);
    hold all
    plot(W_target,EW,'.')
    hold all
    EW2=EW*std(W_target(:))/std(EW(:));
    plot(W_target,EW2,'.')

    legend('x=y','EW','EW2')
    xlabel('True weights')
    ylabel('Estimated weights')
    [R,C,Z,S] = GetWeightsErrors( W_target,EW );
    [R2,C2,Z2,S2] = GetWeightsErrors( W_target,EW2 );

    title({[' EW R =' num2str(R) ', EW2 R =' num2str(R2)]; ...
         [' EW C =' num2str(C) ', EW2 C =' num2str(C2)]; ...
         [' EW Z =' num2str(Z) ', EW2 Z =' num2str(Z2) ];...
         [' EW S =' num2str(S) ', EW2 S =' num2str(S2) ]});
    hold off

pause(1e-6);
    %% Save
    file_name=GetName(params);  %need to make this a meaningful name
    save([file_name(1:end-4) '_SingleNeuron.mat'],'EW','Ebias','W','W_target','bias',...
        'quality','error_rates','Cxx','Cxx_diag_single','Cxy',...
        'Cxy_single','rates','rates_single','params'); 
end
