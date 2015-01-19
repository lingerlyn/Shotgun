function spikes = Calcium2Spikes(Y,noise_method)
% Y - Tx1 or 1xT calcium trace
% noise_method - method used to infer noise power. Can be 'arpfit','psd' or 'mcmc'
% output:
% spikes - de-convolbed spikes

%% pre-process data
addpath('utilities/'); 
if size(Y,2)>1 
    if size(Y,1)==1
        Y=Y'; %makes sure that Y is a column vector
    else
        error('Y must be a vector!')
    end
end
Y =2*Y/max(Y); %re-normalize Y

%% initialize parameters
T = length(Y);
P = arpfit(Y',2);
Gx = @(x,mode) G_mat(x,mode,T,P.g,1);

switch noise_method
    case 'arpfit'
        noise_residual = P.sn*sqrt(T);  
    case 'psd'
        [ ~, noise_residual ] = GetSn(Y);
    case 'mcmc'
        error('mcmc noise method not written yet - add this!')
    otherwise
        error('unknown noise method')
end

%% Re-weighted L1
addpath('utilities/spgl1-1.8'); 
resparse_eps = 1e-6;
N_resparse=3; % number of times to resparse the solution using re-weighted l1-norm

for ii = 1:N_resparse
    if ii>1        
        wt = [1./(spikes+resparse_eps);1];
        options = spgSetParms('project',@NormL1NN_project ,'primal_norm', @NormL1NN_primal,'dual_norm',@NormL1NN_dual,'verbosity',0,'weights',wt);
    else
        options = spgSetParms('project',@NormL1NN_project ,'primal_norm', @NormL1NN_primal,'dual_norm',@NormL1NN_dual,'verbosity',0);
    end    
    tic
    [sp_temp,~,~,~] = spg_bpdn( Gx, Y, noise_residual,options);
    toc
    spikes = sp_temp(1:T);
end

spikes=spikes/max(spikes(:)); %make maximum eqaul 1 so we  can interpet this as spikes


end