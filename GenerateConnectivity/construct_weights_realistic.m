function W = construct_weights_realistic(N,inhib_frac)
%CONSTRUCT_WEIGHTS returns the weight matrix such that
% connectivity is randomly drawn from some diftribution

%INPUT
% N = total number of neurons
% spar = sparsity level
% inhib_frac = fraction of  inhbitory neurons in the population
% OUTPUT
% N x N weight matrix

thresh_max=3;   % maximum synaptic weight
thresh_min=0.01; % minimum synaptic weight
diagonal_weight=2;
W = -diag(diagonal_weight+sqrt(diagonal_weight/10)*randn(N,1)); % diagonal negative
% W = -thresh*diag(rand(N,1)); % diagonal negative

sgn_array=ones(N,1);
sgn_array(1:(1/inhib_frac):end)=-1;

%%  Physiological distance dependent connection probability
D=3;
centers=rand(N,D); %assume neurons are on a D lattice                
dist=zeros(N,N); %distances between neurons
for ii=1:N
    dist(:,ii)=sqrt(sum((mod(bsxfun(@plus,centers(ii,:),-centers)+0.5,1)-0.5).^2,2)); %distance metric - assumes all neurons are on a 1D box with cyclic boundary conditions  
end
rho=1e-4; % neuronal density in 1/(micrometer)^3) from Braitenberg and Schuz, 1991
L=(N/rho).^(1/D); %side length field of view box 
lambda=200; % length constant in micrometers from Perin29032011, Figure 1E
p0=0.2; %initial connectivty probability at d=0 from Perin29032011,Figure 1E 
p=p0*exp(-dist*L/lambda); %connection probability matrix

%% connection strengths
% m_EI=1;
% m_II=0.8*m_EI;
% m_EE=m_EI*inhib_frac/(1-inhib_frac);
% m_IE=0.8*m_EE;

m_EI=1;
m_II=0.8*m_EI;
m_EE=m_EI*inhib_frac/(1-inhib_frac);
m_IE=0.8*m_EE;
v=m_EE.^2;      

% m_EI=5;
% m_II=4;
% m_EE=10^(-0.31);
% m_IE=0.8*10^(-0.31);

for ii=1:N
    sgn=sgn_array(ii);

%     p=GetProb(N,spar,ii);
    conn=find(rand(N,1)<p(:,ii));
    conn(conn==ii)=[]; % since we already have the diagonal
    
    if sgn==1
        conn_type=sgn_array(conn)>0;
        m=m_EE*conn_type+m_IE*(1-conn_type);       
        mu=log(m.^2+(m.^2+v));
        si=sqrt(log(1+v./m.^2));        
        % make sure weights don't exceed cutoff        
        ind=ones(length(conn),1)>0.5;
        while any(ind)                    
            W(conn(ind),ii)=sgn*exp(mu(ind)+randn(length(conn(ind)),1).*si(ind));
            ind=or((W(conn,ii)>thresh_max),(W(conn,ii)<thresh_min));
        end
    elseif sgn==-1
        conn_type=sgn_array(conn)>0;
        m=m_EI*conn_type+m_II*(1-conn_type);
        ind=ones(length(conn),1)>0.5;
        while any(ind)                    
            W(conn(ind),ii)=sgn*m(ind)+randn(length(conn(ind)),1).*sqrt(v);
            ind=or((sgn*W(conn,ii)>thresh_max),(sgn*W(conn,ii)<thresh_min));
        end
        
%         conn_type=sgn_array(conn)<0;
%         m=m_EI*(1-conn_type)+m_II*conn_type;
%         v=m_EE;          
%         mu=log(m.^2+(m.^2+v));
%         si=sqrt(log(1+v./m.^2));        
%         % make sure weights don't exceed cutoff        
%         ind=ones(length(conn),1)>0.5;
%         while any(ind)                    
%             W(conn(ind),ii)=sgn*exp(mu(ind)+randn(length(conn(ind)),1).*si(ind));
%             ind=or((-W(conn,ii)>thresh_max),(-W(conn,ii)<thresh_min));
%         end
        
    end
end

%% Sort W according to distances
chosen_neurons=zeros(N,1)>1; 
current_neuron=1; %first chosen neuron
sorted_neurons=[];
for ii=1:N
    sorted_neurons(end+1)=current_neuron; %#ok
    chosen_neurons(current_neuron)=1>0;
    d=dist(current_neuron,:);
    d(chosen_neurons)=inf;
    [~,current_neuron]=min(d);    
end
W=W(sorted_neurons,sorted_neurons);

%% put all the inhibitory neurons in the end
% W_temp=W;
% W_temp(eye(N)>0.5)=0;
% inhib=find(sum(W_temp,1)<0);
% excit=find(sum(W_temp,1)>=0);
% sorted_neurons=[excit inhib];
% W=W(sorted_neurons,sorted_neurons);

end



