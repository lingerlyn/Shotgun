function c=SetBiases(A,target_rates)
% Estimate the biases needed for each neuron to fire at its target rate.

N=length(A);

zN=sparse(N,N);
W=[zN .5*A zN; .5*A' zN .5*A; zN .5*A' zN];
W0=W(N+1:end,N+1:end);
WT=W(1:2*N,1:2*N);

%some index vectors
nn=(1:N)'; nn2=((N+1):(2*N))'; nn3=((2*N+1):(3*N))';
c=zeros(N,1);
for k=1:N %set the offset term c to produce the desired firing rate
    ind_k=[1:k-1 k+1:N]';
    logit_k=log(target_rates(k)/(1-target_rates(k)));
    c(k)=(logit_k-2*(W(N+k,[nn;N+ind_k;nn3])*[target_rates;target_rates(ind_k);target_rates]));
end;
