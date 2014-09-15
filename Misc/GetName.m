function file_name = GetName( params )
% generate file name for saves
%   Detailed explanation goes here

T=params.spike_gen.T;
N=params.connectivity.N;
obs=params.spike_gen.sample_ratio;
sample_type=[];
if strcmp(params.spike_gen.sample_type,'fixed_subset')
    sample_type='_fixed_subset';
end


file_name=fullfile('Results',['Run_N=' num2str(N) '_obs=' num2str(obs) '_T=' num2str(T) sample_type '.mat']);

end

