function spikes = GetSpikes(W,bias,T,seed,type)
% This function simulates a network with parameters
% W - network connectivity (NxN)
% bias - bias (Nx1)
% T - simulation duration (scalar)
% seed - random seed
% and outputs 
% spikes - network activity (NxT)

switch type
    case 'linear'
        spikes=network_simulation_linear(W,bias,T,seed);
  case 'sign'
        spikes=network_simulation_sign(W,bias,T,seed);
    case 'logistic'
        spikes=network_simulation_logistic(W,bias,T,seed);
end

end

