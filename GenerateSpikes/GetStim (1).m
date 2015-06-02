function stim=GetStim(N_stim,T,stim_type)
%GETSTIM Summary of this function goes here
%   Detailed explanation goes here
stim=zeros(N_stim,T);
if N_stim==0    
    return
end

switch stim_type
    case 'none'
        disp('stimulus=0');
    case 'pulses'
        T_pulse=1e3;
        duty_cycle=0.5;
        pulse_mag=1;
        pulse_std=0.1*pulse_mag;
        
        t_cycle=mod(1:T,T_pulse)>T_pulse*duty_cycle;
        stim=bsxfun(@times,(pulse_mag)+pulse_std*randn(N_stim,T),t_cycle);     
     case 'delayed_pulses'
        pulse_mag=1;
        
        for ii=1:N_stim
            stim(ii,:)=pulse_mag*((mod(1:T,N_stim)+1)==ii);
        end
     case 'Markov'
        state_num=N_stim;
        pulse_mag=30*(rand(state_num,1)-0.5);
        pulse_std=15;
        A=0.00004*diag(rand(state_num,1)); %rate of each state
        
        state=zeros(state_num,1);
        state(randi(state_num))=1;        
        
        transition_mat=A*rand(state_num);
        transition_mat(eye(state_num)>0.5)=0;
        transition_mat(eye(state_num)>0.5)=1-sum(transition_mat,1);
        
        for tt=1:T
            p_next=transition_mat*state;            
            state=diff([0; cumsum(p_next)>rand]);
            stim(:,tt)=(pulse_mag+pulse_std*randn).*state;
        end
     case 'sine'
        T_pulse=1e1;
        pulse_mag=0.5;
        pulse_std=0.5;
        
        stim=bsxfun(@times,(pulse_mag+pulse_std*randn(N_stim,T)),sin(2*pi*(1:T)/T_pulse));
    otherwise
        error('unknown stim type');
end
end

