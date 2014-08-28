function spikes = GetSpikes(W,bias,T,T0,seed,type)
% This function simulates a network with parameters
% W - network connectivity (NxN)
% bias - bias (Nx1)
% T - simulation duration (scalar)
% T0 - burn-in time (scalar) - time to wait so network activity becomes stationary
% seed - random seed
% and outputs 
% spikes - network activity (NxT)

switch type
    case 'linear'
        spikes=network_simulation_linear(W,bias,T,T0,seed);
    case 'linear_reg'
        spikes=regression_simulation_linear(W,bias,T,seed);
    case 'sign'
        spikes=network_simulation_sign(W,bias,T,T0,seed);
    case 'Poisson'
        spikes=network_simulation_Poisson(W,bias,T,T0,seed);
    case 'logistic'
        spikes=network_simulation_logistic(W,bias,T,T0,seed);
end

if any(~isfinite(spikes(:)))
    error('spikes contain non-finite or not defined values')
end

end

