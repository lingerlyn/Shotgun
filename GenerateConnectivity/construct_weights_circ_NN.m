function W=construct_weights_circ_NN(N,NN_range)
% constuct a circular network
% N - number of neurons
% NN_range - distance to connected nearset neighbors
    
    Magnitude=5;

       
    if NN_range>N
        error('NN_range>N !!')
    end
    
    W=zeros(N,N);
    sigma=NN_range/2;
    
    for ii=1:N
        if (ii-NN_range<1)||(ii+NN_range>N)
            NN_indices=[(mod(ii-1-NN_range,N)+1):N,1:(mod(ii-1+NN_range,N)+1)];                    
        else
            NN_indices=(mod(ii-1-NN_range,N)+1):(mod(ii-1+NN_range,N)+1);
        end
        W(ii,NN_indices)=Magnitude*normpdf(-NN_range:NN_range,0,sigma)*rand;
%         W(ii,NN_indices)=Magnitude*(rand(1,2*NN_range+1))/sqrt(2*NN_range+1);
%         W(ii,ii) = -Magnitude/sqrt(2*NN_range+1); % take diagonal weights to cancell all the rest
    end
    
    W(1:(N+1):end) = -Magnitude*ones(N,1)*rand/sqrt(2*pi)/sigma; % take diagonal weights to be negative 
    
    

end