function file_name = GetName( params )
% generate file name for saves
%   Detailed explanation goes here

T=params.spike_gen.T;
N=params.connectivity.N;
obs=params.spike_gen.sample_ratio;
sample_type=[];
if strcmp(params.spike_gen.sample_type,'fixed_subset')
    sample_type='_fixed_subset';
elseif strcmp(params.spike_gen.sample_type,'continuous')
    sample_type='_continuous';
end
    
if params.spike_gen.N_stim>0
    stim_str=['_N_stim=' num2str(params.spike_gen.N_stim)];
else
    stim_str=[];
end

if strcmp(params.spike_gen.neuron_type,'LIF')   
    neuron_str=['_LIF'];
else
    neuron_str=[];
end

if strcmp(params.conn_est_flags.est_type,'Gibbs')   
    est_type_str=['_Gibbs'];
elseif strcmp(params.conn_est_flags.est_type,'FullyObservedGLM')   
    est_type_str=['_FullyObservedGLM'];
elseif strcmp(params.conn_est_flags.est_type,'Cavity')   
    est_type_str=['_Cavity'];
else
    est_type_str=[];
end


file_name=fullfile('Results',['Run_N=' num2str(N) '_obs=' num2str(obs) '_T=' num2str(T) stim_str sample_type neuron_str est_type_str '.mat']);

end

