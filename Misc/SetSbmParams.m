function sbm=SetSbmParams(N,weight_scale)

    DistDep=1; %distance dependent connectivity
    Realistic=1; %Adhere to Dale's law and add a negative diagonal
    str_var=.005; %variance of block weights
    blockFracs=[1/2;1/2];
    
if DistDep        
    sigm=@(z) 1./(1+exp(-z));
    distfun=@(a,b,x) sigm(a*x+b);
%         neuron_positions=linspace(.5*1/N,1-.5*1/N,N)';
    neuron_positions=linspace(0,1,N);
    neuron_positions=neuron_positions(randperm(N)); %shuffle neuron positions
    distdep_a=-15;
    distdep_b=1.3;
    pconn=[];

else
    pconn=spar*ones(length(blockFracs));
    distfun=[];
    neuron_positions=[];
    distdep_a=[];
    distdep_b=[];
end

if Realistic
    self_inhibition=-1*weight_scale; %weight of the diagonal
    idenTol=.25; %Tolerance for classifying neurons as exc or inh
else
    self_inhibition=[];
    idenTol=[];
end

% Set up the SBM mean matrix

    nblocks=length(blockFracs);
    block_means=ones(nblocks).*[ones(nblocks,1) -2*ones(nblocks,1)];

    %increase the rank of MeanMatrix by scaling the values a little
    block_means(1,1)=block_means(1,1)*1.1;
    block_means(1,2)=block_means(1,2)*.9;

    MeanMatrix=GetBlockMeans(N,blockFracs,block_means)*weight_scale;
    if Realistic
        MeanMatrix(~~eye(N))=self_inhibition; 
    end
    
    if DistDep
        neuron_distances=abs(repmat(neuron_positions,[1,N])-repmat(neuron_positions',[N,1]));
        abs(mod(bsxfun(@plus,neuron_positions',-neuron_positions)+0.5,1)-0.5);
        fd=distfun(distdep_a,distdep_b,neuron_distances); %distance-dependent connection probs
    else
        neuron_distances=[];
        fd=[];
    end
    
    sbm=v2struct(Realistic,DistDep,blockFracs,nblocks,str_var,block_means,MeanMatrix,self_inhibition,neuron_positions,...
        neuron_distances,fd,distfun,distdep_a,distdep_b,pconn,idenTol);


end