function [ Y,spikes] = Spikes2Calcium2Spikes( true_spikes )
%SPIKES2CALCIUM2SPIKES Summary of this function goes here
%   Detailed explanation goes here

% observation model
sn=0.2;
amp=1;
b=0; 
g=0.95;
noise_mode=0;
Y = Spikes2Calcium(true_spikes,g,b,amp,sn,noise_mode);

order=1;
% P=GetParams(Y,order,'keep','arpfit');
P= GetParams(Y,order,'psd','arpfit');

[spikes,b] = Calcium2Spikes(Y,P);
thresh=max(spikes,[],2)/2;
spikes=double(bsxfun(@gt,spikes,thresh));

end

