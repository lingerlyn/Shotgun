function [p,dist]=GetDistProb(centers)
    N=size(centers,1);
    D=size(centers,2);
    dist=zeros(N,N); %distances between neurons
    for ii=1:N
        dist(:,ii)=sqrt(sum((mod(bsxfun(@plus,centers(ii,:),-centers)+0.5,1)-0.5).^2,2)); %distance metric - assumes all neurons are on a 1D box with cyclic boundary conditions  
    end
    rho=1e-4; % neuronal density in 1/(micrometer)^3) from Braitenberg and Schuz, 1991
    L=(N/rho).^(1/D); %side length field of view box 
    lambda=200; % length constant in micrometers from Perin29032011, Figure 1E
    p0=0.2; %initial connectivty probability at d=0 from Perin29032011,Figure 1E 
    p=p0*exp(-dist*L/lambda); %connection probability matrix
end