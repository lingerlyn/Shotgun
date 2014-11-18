function observations=SampleSpikes(N,T,sampled_ratio,sample_type,N_stim,seed)
% inputs:
% spikes - NxT spikes matrix
% sample_ratio - perscent of observed neurons
%  sample_type - see below

% outputs:
% ind - indices of observed spikes


unsampled_ratio=1-sampled_ratio;
ind=zeros(N,T)>1;

stream = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(stream);

    switch sample_type
        case 'spatially_random'
            for tt=1:T
                temp=randperm(N);
                ind(temp(1:round(N*unsampled_ratio)),tt)=1>0; %this randomization insures a at least sample_ratio observed at each time
            end 
        case 'prob'
            addpath('Misc') % so we can use the GetProb function            
            spar=0.2; %for now I hardwire this. later, make this one of the function's inputs
            temp=randperm(N);
            ind(temp(1:round(N*sampled_ratio)),1)=1>0; %note we use here sampled_ration and not unsampled_ratio
            for tt=2:T
                ind_num=find(ind(:,tt-1));
                p=GetProb(N,spar,ind_num);
                p=bsxfun(@times,p,1./sum(p,1));                
                next_ii=SampleProb(p);
                % if next indices have duplicates, replace these duplicates with random choces
                next_ii=sort(next_ii);                
                duplicates=[diff(next_ii)==0 0]>0.5;     
                if ~isempty(duplicates)
                    unique_next_ii=unique(next_ii);
                    Complmentary_set=1:N;
                    for jj=1:length(unique_next_ii)
                        Complmentary_set(Complmentary_set==unique_next_ii(jj))=0;
                    end
                    Complmentary_set(Complmentary_set==0)=[];
                    temp=Complmentary_set(randperm(length(Complmentary_set)));
                    next_ii(duplicates)=temp(1:sum(duplicates));         
                end
                ind(next_ii,tt)=1>0; %this randomization insures a at least sample_ratio observed at each time
                
            end 
            ind=~ind; %we will inverse this back again in the end
        case 'continuous'
            for tt=1:T                
                ind(1+mod(tt+(1:round(N*unsampled_ratio)),N),tt)=1>0;
            end      
        case 'fully_random'
            ind=rand(N,T)<unsampled_ratio;
        case 'fixed_subset'
            ind((ceil(N*sampled_ratio)+1):end,:)=1;
        case 'random_fixed_subset'
            temp=randperm(N);
            ind(temp(ceil(N*sampled_ratio):end),:)=1;
        case 'random_NN_blocks'
            K=floor(unsampled_ratio*N/2);
            for tt=1:T                
                temp=mod((-K:K)+randi(N),N)+1;
                ind(temp,tt)=1>0; %this randomization insures a at least sample_ratio observed at each time
            end
        case 'serial_scan'
            block_size=floor(unsampled_ratio*N);
            block_num=ceil(1/unsampled_ratio);
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
    
    observations=[~ind; ones(N_stim,T)];
end

function x=SampleProb(p)
% draw a single sample x from a discrete probability distribution p (column vector)
% if p is a matrix - do this for each column vector

r = rand(1,size(p,2));
x=sum(bsxfun(@gt, r,cumsum(p)))+1;

end