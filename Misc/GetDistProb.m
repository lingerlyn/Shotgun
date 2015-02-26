function [p,dist]=GetDistProb(centers,spar)
    % input:
    % centers - locations in [0,1]^D box
    % spar - sparsity level. if empty, we use default 0.2
    % output:
    % dist - distances between neurons, in micrometes
    % p - probabilites of connections

    N=size(centers,1); %number of neurons
    D=size(centers,2); %number of dimensions
    rho=1e-4; % neuronal density in 1/(micrometer)^3) from Braitenberg and Schuz, 1991
    L=(N/rho).^(1/D); %side length field of view box 
    lambda=200; % length constant in micrometers from Perin29032011, Figure 1E
    p0=0.2; %initial connectivty probability at d=0 from Perin29032011,Figure 1E 
    
    if ~isempty(spar)
        p0=spar;
    end
    
    dist=zeros(N,N); %distances between neurons    
    for ii=1:N
        dist(:,ii)=L*sqrt(sum((mod(bsxfun(@plus,centers(ii,:),-centers)+0.5,1)-0.5).^2,2)); %distance metric - assumes all neurons are in a box with cyclic boundary conditions  
    end
    p=p0*exp(-dist/lambda); %connection probability matrix
end