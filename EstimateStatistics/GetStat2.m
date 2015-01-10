function [ CXX_filtered, CXY_filtered,W,rates_filtered,obs_count] = GetStat2( sampled_spikes,observations,filter_list,glasso,restricted_penalty,pos_def,sparsity,true_W)
% inputs:
% sampled_spikes -NxT observed spikes (with zeros in unobsereved samples)
% observations - NxT matrix of which timebins were observed
% filter_list - 2xNxL cell array with lists of poles and zeros for filters bank
% glasso - flag that indicate whether or not use glasso
% restricted_penalty - flag that indicates whether or not to use a restricted l1 penality in lasso (only on parts of the inv_COV matrix)
% pos_def - restict Covariance matrix to be positive semi definite
% sparsity-  the average nnz in the inv(CXX) matrix.  if empty, test using true matrix
% true_W - for  testing purposes

% outputs:
% Cxx_filtered  - NxNxL Array of estimated filtered spike covariance at a single timestep, one for each of the L filters
% Cxy_filtered - NxNxL matrix of estimated cross-covariance of filtered spikes to spikes between two adjacent timesteps, one for each of the L filters
% W - infered weights matrix
% rates_filtered - NxL mean firing rates, after filtering with each filter
% obs_count - observation count for each neuron


% internal parameters
Tol=1e-6; %numerical tolerance for glasso
msg=1; %verbosity of glasso algorithm
maxIter=100;  %glasso max iterations
Tol_sparse=1e-2; % tolerance for sparsity of W
thresh=inf; %max size of spikes array for which we use effiecient computation of sufficient statistics

[N, T] = size(sampled_spikes);
L=size(filter_list,3);
filtered_spikes=zeros(N,T,L);
filter_square=zeros(N,N,L);
filter_sum=zeros(N,L);
delta_function_array=zeros(T,1); %does not have to be T  long - just longer then the filter
delta_function_array(2)=1;

for ll=1:L
    for nn=1:N
        a=filter_list{1,nn,ll}; %polynom of transfer function denomenator    
        b=filter_list{2,nn,ll}; %polynom of transfer function numerator
        filtered_spikes(nn,:,ll)=filter(b,a,sampled_spikes(nn,:));
        temp1=filter(b,a,delta_function_array);                
        filter_sum(nn,ll)=sum(temp1);
        for kk=1:N                
                c=filter_list{1,kk,ll}; %polynom of transfer function numerator
                d=filter_list{2,kk,ll}; %polynom of transfer function denomenator    
                temp2=filter(d,c,delta_function_array);
                filter_square(nn,kk,ll)=sum(temp1.*temp2);                
        end
    end
end

mYn=sum(observations,2);
mY=sum(sampled_spikes,2);
rates=mY./(mYn+eps); %estimate the mean firing rates
rates_filtered=zeros(N,L);
for ll=1:L
    rates_filtered(:,ll)=filter_sum(:,ll).*rates;
end

XX_filtered=zeros(N,N,L);
XY_filtered=zeros(N,N,L);
CXX_filtered=zeros(N,N,L);
CXY_filtered=zeros(N,N,L);

%initialize the sufficient statistics arrays
observations=double(observations);
if N*T<thresh
        XX=sampled_spikes*sampled_spikes'/T;
        XXn=observations*observations'/T;
        XYn=observations(:,1:(end-1))*(observations(:,2:end))'/T; %or mYn*mYn'/T
    for ll=1:L
        XX_filtered(:,:,ll)=filtered_spikes(:,:,ll)*filtered_spikes(:,:,ll)'/T;
        XY_filtered(:,:,ll)=filtered_spikes(:,1:(end-1),ll)*(sampled_spikes(:,2:end))'/T;        
    end
else %haven't written this part yet - after writing this, change thresh to 1e8
    XX=zeros(N);
    XXn=zeros(N);
    XY=zeros(N);
    XYn=zeros(N);
    g = find(observations(:,1));
    sg = double(sampled_spikes(g,1));
    disp(['Gathering sufficient statistics...'])
    for tt = 2:T
        f = g;
        sf = sg;
        g = find(observations(:,tt));
        sg = double(sampled_spikes(g,tt));

        XX(f,f)=XX(f,f)+sf*sf';
        XXn(f,f)=XXn(f,f)+1;

        XY(f,g)=XY(f,g)+sf*sg';
        XYn(f,g)=XYn(f,g)+1;
        if ~mod(tt,T/10)
            disp([num2str(100*tt/T,2) '%'])
        end
    end
end


%% Calculate sufficient stats

for ll=1:L
    CXX_filtered(:,:,ll)=(XX_filtered(:,:,ll)-filter_square(:,:,ll).*XX(:,:,ll).*(1-XYn(:,:,ll)./(XXn(:,:,ll)+eps)))./(XYn(:,:,ll)+eps)-rates_filtered(:,ll)*rates_filtered(:,ll)'; %estimate the covariance
    % CXX((XXn<10))=0;%set elements to zero that haven't been observed sufficiently
    CXY_filtered(:,:,ll)=XY_filtered./(XYn+eps)-rates_filtered(:,ll)*rates(:,ll)'; %estimate cross-covariance
    % CXY((XYn<10))=0;
    COV = [CXX_filtered(:,:,ll) CXY_filtered(:,:,ll); CXY_filtered(:,:,ll)' CXX_filtered(:,:,ll)];

    if(any(eig(COV)<0))
        disp('COV is not positive semidefinite;')
        if pos_def %positive semidefinite projection, so we won't have problems with glasso
            disp('correcting...');
            [v,d]=eig(COV);
            X0=v*spdiags(max(diag(d),0),0,2*N,2*N)*v';
            COV=X0;
        end
    end

    inv_COV=COV\eye(2*N);

    %% Use glasso to sparsify stats
    if glasso==1    
        disp('starting glasso...')

        addpath(fullfile('EstimateStatistics','QUIC')) %mex files for QUIC glasso implementation

        lambda_high=1; %maximum bound for lambda
        lambda_low=1e-4;  %minimum bound for lambda
        if restricted_penalty %more accurate, but much slower
            regularization_mat=ones(2*N);
            regularization_mat(1:N,1:N)=0; %remove penalty with upper left block
            regularization_mat(eye(2*N)>0.5)=0;  %remove penalty from diagonal
            regularization_mat([zeros(N), eye(N); zeros(N), zeros(N)]>0.5)=0;  %remove penalty from diagonal of upper right block
            regularization_mat([zeros(N), zeros(N); eye(N), zeros(N)]>0.5)=0;  %remove penalty from diagonal of lower left block         
        else
             regularization_mat=1;
        end


        if isempty(sparsity) % in case we "cheat" choouse optimal lambda
            func= @(x) OptFun( x,COV,regularization_mat, Tol, msg, maxIter ,true_W);
            lambda_opt=fminbnd(func,lambda_low,lambda_high);
            lambda_mat=lambda_opt*regularization_mat;
            [inv_COV_res, COV_res, opt, cputime, iter, dGap] = QUIC('default', COV, lambda_mat, Tol, msg, maxIter);

        else %set lambda according to estimated sparsity level of W, using binary search
            loop_cond=1;  %flag for while llop
            flag_first=1; % flag for first iteration of loop    

            while  loop_cond
                if sparsity==1
                    lambda=0;
                else
                    lambda=(lambda_high+lambda_low)/2;
                end
                    lambda_mat=lambda*regularization_mat;
                    if flag_first
                        [inv_COV_res, COV_res, opt, cputime, iter, dGap] = QUIC('default', COV, lambda_mat, Tol, msg, maxIter);
                    else %use warm start
                        flag_first=0;
                        [inv_COV_res, COV_res, opt, cputime, iter, dGap] = QUIC('default', COV, lambda_mat, Tol, msg, maxIter,inv_COV_res, COV_res);
                    end

                    % W block
                    temp=inv_COV_res((N+1):end,1:N);
                    temp(eye(N)>0.5)=0;

                    sparsity_measure=mean(~~temp(:)); 
                    cond=sparsity_measure<sparsity;
                    loop_cond=(abs(sparsity_measure-sparsity)/sparsity >  Tol_sparse);
                     if sparsity==1
                        loop_cond=0;
                     end

                if cond
                    lambda_high=lambda
                else
                    lambda_low=lambda
                end
            end
        end

            CXY = COV_res(1:N,N+1:end);%%%%%
            CXX = COV_res(1:N,1:N);%%%%%
            W=-inv_COV_res((N+1):end,1:N);            
        else
            CXY = COV(1:N,N+1:end);%%%%%
            CXX = COV(1:N,1:N);%%%%%
            W=-inv_COV((N+1):end,1:N); % W estimate without glasso        
    end
end
obs_count=mYn;
end

