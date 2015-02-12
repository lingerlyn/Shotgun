function samples = SampleDist(p,sample_num)
% inputs:
% p - a probability matrix (sum of columns is 1)
% sample_num - number of samples to draw
% output
% samples - a sample_num size vector of samples drawn from p

% assert(abs(sum(p)-1)<1e-6,'p vector must sum to 1');
assert(all(p(:)>=0),'p vector must no negative');
assert(sample_num>0,'sample_num must be positive');
[N,L]=size(p);
samples=zeros(N,sample_num);

for nn=1:N
    temp=[0 cumsum(p(nn,:))]';
    rand_seeds=rand(1,sample_num);
    samples_binary=diff(bsxfun(@gt, rand_seeds, temp),[],1)<0;
    [indices,~]=find(samples_binary);
    samples(nn,:)=(indices-1)/L; % normalize samples to be in [0 1]
end

end

