function sampled_spikes=SampleSpikes(spikes,sample_ratio,sample_type,seed)
% inputs:
% spikes - NxT spikes matrix
% sample_ratio - perscent of observed neurons
%  sample_type - see below

% outputs:
% ind - indices of observed spikes

N=size(spikes,1);
T=size(spikes,2);

ind=zeros(N,T)>1;

stream = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(stream);

    switch sample_type
        case 'spatially_random'
            for tt=1:T
                temp=randperm(N);
                ind(temp(1:round(N*sample_ratio)),tt)=1>0; %this randomization insures a at least sample_ratio observed at each time
            end            
        case 'fully_random'
            ind=rand(N,T)<sample_ratio;
        case 'fixed_subset'
            ind(1:N*sample_ratio,:)=1;
        case 'random_NN_blocks'
            K=floor(sample_ratio*N/2);
            for tt=1:T                
                temp=mod((-K:K)+randi(N),N)+1;
                ind(temp,tt)=1>0; %this randomization insures a at least sample_ratio observed at each time
            end
        case 'serial_scan'
            block_size=floor(sample_ratio*N);
            block_num=ceil(1/sample_ratio);
            for tt=1:T
                if block_num>2
                    block_start=mod(tt,block_num)*block_size;
                else
                    block_start=mod(tt,block_num)*(N-block_size);
                end
                temp=block_start+(1:block_size);
                ind(temp,tt)=1>0; %this randomization insures a at least sample_ratio observed at each time
            end
        otherwise
            error('unknown sample_type!!');
    end
    
    sampled_spikes = spikes;
    sampled_spikes(ind)=nan;
end