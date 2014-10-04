function EW=OMPNL(Cxx,Cxy,spar,rates)%,M,lambda)
%add penalties later.
N=length(Cxx);
h=-rates.*log(rates)-(1-rates).*log(1-rates); %entropy

maxCGiter=100;

EW=zeros(N);
params.Cxx=Cxx;
params.Cxy=Cxy;
params.h=h;

for i=1:N
    w=zeros(N,1);
    S=[]; %support
    Sc=(1:N)'; %complement of support
    params.supp=[];
    params.Cxy=Cxy(:,i);
    params.h=h(i);
    params.w=w;
    
    while mean(~~w)<spar
 
    %find optimal supports and their values given current w

    nSc=length(Sc);
%     idx_val=zeros(nSc,1);
    fX=zeros(nSc,1);

    for j=1:nSc
        params.idx=Sc(j); %index of w_i
    
%         idx_val(j)=minimize(0,'GetLLandGradGivenOthers',maxCGiter,params);
        [~,ff]=minimize(0,'GetLLandGradGivenOthers',maxCGiter,params);
        fX(j)=ff(end);
    end
    %find argmax
    fX=-fX; %since we minimized its negation.
    [~,idx]=max(fX);

    %update support
    S=sort([S;Sc(idx)]);
    Sc(idx)=[];
    params.supp=S;
    
    %Get best w given the new support
%     ww=minimize(w(S),'GetLLandGrad',maxCGiter,params); %use previous w as guess? might be better to reset to zeros.
    ww=minimize(zeros(length(S),1),'GetLLandGrad',maxCGiter,params);

    w(S)=ww;
    params.w=w; %update current estimate of w.
    end
    EW(i,:)=w;
end


end


