function Y = Spikes2Calcium(spikes,g,b,amp,noise_std)
%GENERATE_Y Generates calcium observation 
% INPUT
% spikes =  N  X T spikes array
% obs_neurons = N x T observation matrix (0-1)
% g = gamma (calcium rates)
% b = baseline
% amp = amplitude
% noise_std = noise std

% OUTPUT:
% Y :  N x T array for the observed calcium for neuron n at time t
addpath('utilities/'); 


[N, T] = size(spikes);
Gx = @(x,mode) G_mat(x,mode,T,g,1);
Y=zeros(N,T);
for ii=1:N
    Y(ii,:) = Gx(spikes(ii,:),1) + amp*noise_std*randn(1,T) + b; 
end
end