function [ Y,spikes] = Spikes2Calcium2Spikes( spikes )
%SPIKES2CALCIUM2SPIKES Summary of this function goes here
%   Detailed explanation goes here

% observation model
sn=0.2;
amp=1;
b=1; 
g=0.95;
noise_mode=0;
Y = Spikes2Calcium(spikes,g,b,amp,sn,noise_mode);

order=1;
% P=GetParams(Y,order,'keep','arpfit');
P= GetParams(Y,order,'psd','arpfit');

[spikes,b] = Calcium2Spikes(Y,P);
spikes=double(spikes>0.5);

end

