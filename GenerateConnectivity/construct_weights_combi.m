function [ W ] = construct_weights_combi(N)
%CONSTRUCT_WEIGHTS Summary of this function goes here
%   Detailed explanation goes here

if N < 20
     error('Take N>=20');
end

total_conn = ceil(0.15*N);
W = diag(unifrnd(-2,0,N,1)); % take diagonal wieghts to be negative 
not_obs_subset = (ceil(N/3)+1):N;
for i = not_obs_subset
    nearby_conn = ceil(total_conn/2);
    for j = 1:nearby_conn
        if i+j<=N
            W(i,i+j)=unifrnd(-2,2,1,1);
        end
    end
    left_conn = total_conn - nearby_conn;
    other_conn = randsample(find(~W(i,:)),left_conn);
    W(i,other_conn)=unifrnd(-2,2,length(other_conn),1);
end

for i = 1:(ceil(N/3))
    conn_unobs_neu = randsample((ceil(N/3)+1):N,2);
    sign = rand<0.5;
    if sign ==0
        sign = -1;
    end
    W(i,conn_unobs_neu(1)) = sign*2;
    W(conn_unobs_neu(1),conn_unobs_neu(2)) = 2;
    W(i+1,conn_unobs_neu(2)) = sign*2;
end