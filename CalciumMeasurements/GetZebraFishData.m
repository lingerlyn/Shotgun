function [data, cent] = GetZebraFishData(subset,N)
% inputs:
% subset - restrict data to certain region? can be 'all','front' or 'back'
% N - number of neurons to take
% output:
% data - (N+1)xT data of extracted calcium activity - N neurons and 1 component of the mean of all the rest
% cent - Nx3 3D locations of neurons

load('C:\Users\Daniel\Copy\Columbia\Research\Group Lasso\Code\PostProcessing\PostProessed_Misha_Data_20140418\Traces&centers_merged_sorted.mat','activity_merged','cent_merged');
t_start=300; %start time - after non-stationarity decline in the beginning of experiment
data=cell2mat(activity_merged);
% data = zscore(data')';
data(1:t_start,:)=[]; %remove tranisent part

if strcmp(subset,'all');
%     x_min=0;
%     x_max=inf;
    ind=1:size(data,2);
elseif strcmp(subset,'front');
    x_min=0;
    x_max=500;    
    ind=find((cent_merged(1,:)>x_min)&(cent_merged(1,:)<x_max));
elseif strcmp(subset,'back');
    x_min=500;
    x_max=inf;
    ind=find((cent_merged(1,:)>x_min)&(cent_merged(1,:)<x_max));
else 
    error('unknown subset type');
end

N=min(N,length(ind)); 
ind=ind(1:N);
all_the_rest=1:size(data,2);
all_the_rest(ind)=[];
external_input=mean(data(:,all_the_rest),2);
external_input=external_input/mean(external_input);

data=[data(:,ind) external_input]';
cent=cent_merged(:,ind)';

end