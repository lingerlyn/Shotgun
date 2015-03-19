    spikes=zeros(N,T);
for nn=1:N
    tic
    Gx = @(x,mode) G_mat(x,mode,T,P{nn}.g,1,P{nn}.z);
    a=ones(T,1);
    m=mean(Gx(a,1));
    num_spikes=round(T*mean(Y(nn,:),2)/m); %later modify this with exponential search
    
    GY=Gx(Y(nn,:),1);
    for ss=1:num_spikes
        [~,ind]=max(GY);    
        spikes(nn,ind)=1;
        temp=zeros(1,T);
        temp(ind)=1;
        GY=GY-Gx(temp,1);
        if ~mod(ss,floor(T/10))
            disp([num2str(100*ss/T,2) '%'])
        end
    end
    toc
end