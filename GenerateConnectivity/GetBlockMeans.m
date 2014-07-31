function str_means_mat = GetBlockMeans(N,blockFracs,str_mean)
% Little helper function that finds the means of an sbm based on the
% 'negative diagonal' pattern.

    labels=zeros(N,1);
    labels( 1:blockFracs(1)*N )=1;
    nTypes=length(blockFracs);
    for k=2:nTypes
        labels( sum(blockFracs(1:k-1))*N+1: sum(blockFracs(1:k-1))*N+blockFracs(k)*N) = k;
    end
    
    str_means_mat=zeros(N);
    for n=1:N
        for nn=1:N
            str_means_mat(n,nn)=str_mean(labels(n),labels(nn));
        end
    end


end