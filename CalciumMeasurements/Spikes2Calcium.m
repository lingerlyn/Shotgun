function Y = Spikes2Calcium(spikes,g,b,amp,noise_std,noise_mode)
%GENERATE_Y Generates calcium observation 
% INPUT
% spikes =  N  X T spikes array
% obs_neurons = N x T observation matrix (0-1)
% g = gamma (calcium rates)
% b = baseline
% amp = amplitude
% noise_std = noise stdev
% noise_mode - if 1, noise is also filtered, if 0, noise is not filtered

% OUTPUT:
% Y :  N x T array for the observed calcium for neuron n at time t
% addpath('utilities/'); 

[N, T] = size(spikes);
Gx = @(x,mode) G_mat(x,mode,T,g,1);
Y=zeros(N,T);
for ii=1:N
    if noise_mode==0
        Y(ii,:) = amp*Gx(full(spikes(ii,:)),1) + noise_std*randn(1,T) + b; 
    elseif noise_mode==1
        Y(ii,:) = amp*Gx(full(spikes(ii,:)) + noise_std*randn(1,T)/2,1)+ noise_std*randn(1,T)/2 + b;         
    end
end
end