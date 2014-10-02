function stim=GetStim(N_stim,T,stim_type)
%GETSTIM Summary of this function goes here
%   Detailed explanation goes here
stim=zeros(N_stim,T);
switch stim_type
    case 'none'
        disp('stimulus=0');
    case 'pulses'
        T_pulse=1e3;
        duty_cycle=0.5;
        pulse_mag=0.7;
        pulse_std=0.1;
        
        t_cycle=mod(1:T,T_pulse)>T_pulse*duty_cycle;
        stim=bsxfun(@times,(pulse_mag+pulse_std*randn(N_stim,T)),t_cycle);     
     case 'delayed_pulses'
        pulse_mag=1;
        
        for ii=1:N_stim
            stim(ii,:)=pulse_mag*((mod(1:T,N_stim)+1)==ii);
        end
     case 'sine'
        T_pulse=1e4;
        pulse_mag=0.5;
        pulse_std=0.5;
        
        stim=bsxfun(@times,(pulse_mag+pulse_std*randn(N_stim,T)),sin(2*pi*(1:T)/T_pulse));
    otherwise
        error('unknown stim type');
end
end

