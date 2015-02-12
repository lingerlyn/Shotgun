function [ Y,spike_output,dt ] = GetTimData(dep)
%GETTIMDATA Summary of this function goes here
%   Detailed explanation goes here
    load('C:\Users\Daniel\Copy\Columbia\Research\Data\Tim_Machado\antidromic-GCaMP6S.mat','antidromicData');
   
    dt = mean(diff(antidromicData{1}.frameOnsets{1}));
    
    Y = antidromicData{1}.traces{dep};
    tt = antidromicData{1}.frameOnsets{dep};
    spikes = antidromicData{1}.stimulusOnsets{dep};
    spiketimes = zeros(size(tt));    

    for i = 1:length(spikes)
        ff = find(tt>spikes(i),1,'first');
        spiketimes(ff) = spiketimes(ff) + 1;
    end
    spike_output = repmat(spiketimes,[size(Y,1) 1]);


end

