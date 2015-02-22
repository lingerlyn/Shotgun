function spikes_out = GetSpikes(W,bias,T,T0,seed,type,N_stim,stim_type,timescale,s0,verbos)
% This function simulates a network with parameters
% W - network connectivity (NxN)
% bias - bias (Nx1)
% T - simulation duration (scalar)
% T0 - burn-in time (scalar) - time to wait so network activity becomes stationary
% seed - random seed
% s0 - (Nx1) initial condition. can be empty. Currently only works for  'logistic_with_history'
% verbos - should we give a progress report?
% and outputs 
% spikes - network activity (NxT)

G=W(1:(end-N_stim),(end-N_stim+1):end);
stim=GetStim(N_stim,T,stim_type);
U_ext=G*stim;    
    
bias=bsxfun(@plus,bias,U_ext);
A=W(1:(end-N_stim),1:(end-N_stim));

switch type
    case 'linear'
        spikes=network_simulation_linear(A,bias,T,T0,seed);
    case 'linear_reg'
        spikes=regression_simulation_linear(A,bias,T,seed);
    case 'sign'
        spikes=network_simulation_sign(A,bias,T,T0,seed);
    case 'Poisson'
        spikes=network_simulation_Poisson(A,bias,T,T0,seed);
    case 'logistic'
        spikes=network_simulation_logistic(A,bias,T,T0,seed);
    case 'logistic_with_history'
        spikes=network_simulation_logistic_with_history(A,bias,T,T0,seed,timescale,s0,verbos);
    case 'logistic_with_delays'
        spikes=network_simulation_logistic_with_delays(A,bias,T,T0,seed);
    case 'LIF'
        spikes=network_simulation_LIF(A,bias,T,T0,seed,timescale,s0,verbos);
    otherwise
            error('unknown neuron type');
end

spikes_out=[spikes; stim];
if mean(spikes(:))<0.1
    spikes_out=sparse(spikes_out);
end

if any(~isfinite(spikes(:)))
    error('spikes contain non-finite or not defined values')
end

end

