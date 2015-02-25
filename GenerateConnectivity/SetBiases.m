function Eb_out=SetBiases(W,target_rates,spike_gen)
% Estimate the biases needed for each neuron to fire at its target rate.

eta=1; % learning rate
T=1e3; % sample simulation time
rates_tol =0.15;
max_iteration=2000;
inhib_target_rates=1.5*target_rates;
thresh=0.2; %maximal firing rate 20Hz (?)

% obtain modify target rates according to neuronal type - inhibitory neurons fire (twice?) more
W_temp=W;
N=size(W,1);
W_temp(eye(N)>0.5)=0;
ind_excitatory=(sum(W_temp,1)>=0)';
ind_inhibitory=(sum(W_temp,1)<0)';
target_rates_with_types=target_rates.*ind_excitatory+inhib_target_rates.*ind_inhibitory;

noise_dist='normal';
switch noise_dist
    case 'normal' % randomize rates according to a normal distriubtion
%          target_rates_with_types=abs(target_rates_with_types+sqrt(target_rates_with_types/10)*randn);

    case 'lognormal'% randomize rates according to a lognormal distriubtion (Figure 3 in "The log-dynamic brain: how skewed distributions affect network operations")
        m=target_rates_with_types;
        v=m;            
        mu=log(m.^2+(m.^2+v));
        si=sqrt(log(1+v./m.^2));      
        ind=ones(N,1)>0.5;
        while any(ind)                    
            target_rates_with_types(ind)=exp(mu(ind)+randn(sum(ind),1).*si(ind));
            ind=target_rates_with_types(:)>thresh;
        end
end 

[~,T0,~,~,seed_spikes,~,N_stim,stim_type, neuron_type,timescale,~]=v2struct(spike_gen);

logit_target_rates=log(target_rates_with_types./(1-target_rates_with_types));
mean_U=W*target_rates_with_types;    
Cxx=diag(target_rates_with_types.*(1-target_rates_with_types)); %we approximate initialy that the spikes are weakly correlated
var_U=diag(W*Cxx*W');   
Eb=sqrt(1+pi*var_U/8).*logit_target_rates-mean_U;
rates_diff=inf;
iteration=1;

verbos=0;
s0=[];
Eb_out=0;

while rates_diff>rates_tol 
    iteration=iteration+1;
    spikes=GetSpikes(W,Eb,T,T0,seed_spikes,neuron_type,N_stim,stim_type,timescale,s0,verbos);
    rates=mean(spikes,2);
%     Eb=Eb-(eta/iteration)*(rates-target_rates_with_types);      
    Eb=Eb-eta*(rates-target_rates_with_types);  
    Eb_out=Eb_out*(1-1/iteration)+Eb/iteration;
    spikes=GetSpikes(W,Eb_out,T,T0,seed_spikes,neuron_type,N_stim,stim_type,timescale,s0,verbos);
    rates_out=mean(spikes,2);
    rates_diff=mean(abs(rates_out-target_rates_with_types))/mean(target_rates_with_types);

    excitatory_rate=mean(mean(spikes(ind_excitatory,:)))
    inhibitory_rate=mean(mean(spikes(ind_inhibitory,:)))
    
    if iteration>max_iteration
        warning('max iteration reached');
        return
    end
end

end

