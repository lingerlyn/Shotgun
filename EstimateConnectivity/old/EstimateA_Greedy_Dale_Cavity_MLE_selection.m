function EW=EstimateA_Greedy_Dale_Cavity_MLE_selection(Cxx,Cxy,rates,spar,lambda,M,idents)

N=length(Cxx);

h=-rates.*log(rates)-(1-rates).*log(1-rates); %entropy

EW=zeros(N);

for a=1:N
    disp(a)
    w=zeros(N,1);
    
    %initialize support
    S=a;
    Sc=[(1:(a-1))';((a+1):N)'];
    
    
    while (nnz(w)/N)<spar
 
        wSc=zeros(length(Sc),1);
        LSc=zeros(length(Sc),1);
        
        %find optimal ws and evaluate the log-likelihoods
        for bb=1:length(Sc)
            b=Sc(bb);
            
            SS=sort([S; b]);
            idx=find(SS==b);
            
            temp=EstimateA_MLE_cavity_restricted(a,SS,Cxx,Cxy,rates);  
            wSc(bb)=temp(idx);

            %look at resulting likelihoods of solutions found
            ww=w;
            ww(b)=temp(idx);
            wnonzero=ww(SS);
            LSc(bb)=wnonzero'*Cxy(SS,a)-h(a)*sqrt(1+pi/8*wnonzero'*Cxx(SS,SS)*wnonzero)-lambda*(wSc(bb)-M(a,b))^2;
 
        end
        
            
        %now form a penalized objective
        obj=LSc-abs(idents(Sc)).*(10^50*(sign(wSc)~=idents(Sc)));
        [~,idx]=max(obj);

        idx=idx(randi(length(idx))); %if the objs are equal pick a random one

        %update support
        S=sort([S;Sc(idx)]);
        Sc(idx)=[];

        w(S)=EstimateA_MLE_cavity_restricted(a,S,Cxx,Cxy,rates);   
    end
    EW(a,:)=w;
end
    
end