function p=GetProb(N,spar,ii)
%% get probability of connection of neuron ii to other neurons on a 1D ring
% lattice with an exponential decay of connection probability
% N -number of neurons
% spar - sparsity of connections
% chosen conenction 

    D=1; %number of dimensions - though we can only handle D=1 for now
    centers=linspace(0,1,N)'; %assume a 1-D lattice
    V=(pi^(D/2))/gamma(1+D/2);
    dist=abs(mod(bsxfun(@plus,centers(ii)',-centers)+0.5,1)-0.5)/(spar*V)^(1/D); %distance metric - assumes all neurons are on a 1D box with cyclic boundary conditions

    p=exp(-3*dist);
end