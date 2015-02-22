function [ Y,cents ,dt ] = GetManolisData(dataname)
%GETTIMDATA Summary of this function goes here
%   Detailed explanation goes here
if ~or(strcmp(dataname,'New'),strcmp(dataname,'Old'))
    error('dataname must "Old" or "New" ');
end

load(['C:\Users\Daniel\Copy\Columbia\Research\Data\Tolias_22.12.2014\Data' dataname '.mat'],'Data');
   
dt = 1/Data.SamplingRate;
cents=Data.Coordinates;
Y=Data.Traces';    
end

