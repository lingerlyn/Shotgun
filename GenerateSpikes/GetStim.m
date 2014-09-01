function stim=GetStim(N_stim,T,stim_type)
%GETSTIM Summary of this function goes here
%   Detailed explanation goes here

switch stim_type
    case 'none'
        stim=zeros(N_stim,T);
    case 'pulses'
        T_pulse=1e2;
        duty_cycle=0.3;
        pulse_mag=0.2;
        pulse_std=0.2;
        
        t_cycle=mod(1:T,T_pulse)>T_pulse*duty_cycle;
        stim=bsxfun(@times,(pulse_mag+pulse_std*randn(N_stim,T)),t_cycle);      
    otherwise
        error('unknown stim type');
end
end

