function X=weighted_semidef_project(X0,D,W,n_iters,p);

%compute argmin_X ||X-D||_W, X \in S
%D and W must be symmetric
%||.||_W is the squared weighted frobenius norm (W is pos weight matrix)
%p is the admm param
%X0 is initial guess

%clear;

test_flag=0;
if(test_flag) %for testing
    tol=1e-4;
    n_iters=20;
    %D=randn(10); D=D+D';
    %W=.01*tol*ones(size(D));
    %V=rand(size(D))<.4; V=V+V'; W(find(V))=10;
    %D(find(~V))=0;
    
    X0=D;
    p=1;
end;

%initialize
Xold=X0;
Zold=X0;
Yold=zeros(size(D));
n=size(D,1);

err=1;
%while err>tol
for iter=1:n_iters
    Xnew=(W.*D-Yold+p*Zold)./(W+p);
    %Znew=argmin_Z sd(Z) - <Z,Y{k}> + .5*p*||Z-Xnew||
    temp=Yold/p+Xnew;
    [v,d]=eig(temp);
    Znew=v*spdiags(max(diag(d),0),0,n,n)*v';
    Ynew=Yold+p*(Xnew-Znew);
    
    err=norm([Xold(:)-Xnew(:);Yold(:)-Ynew(:);Zold(:)-Znew(:)]);
    Xold=Xnew;
    Yold=Ynew;
    Zold=Znew;
    
    disp(sum(W(:).*(Xnew(:)-D(:)).^2))
    
    if(test_flag)
        figure(1);
        minc=min([Xnew(:);D(:);Znew(:)]);
        maxc=max([Xnew(:);D(:);Znew(:)]);
        mmc=[minc maxc];
        
        subplot(231); imagesc(Xnew); caxis(mmc); colorbar;
        subplot(232); imagesc(Znew); caxis(mmc); colorbar; title(err)
        subplot(233); imagesc(Ynew); colorbar;
        subplot(234); imagesc(D); caxis(mmc); colorbar;
        subplot(235); imagesc(W); colorbar;
        pause(.01);
    end;
end;
[v,d]=eig(Xnew);
X=v*spdiags(max(diag(d),0),0,n,n)*v';

disp(sum(W(:).*(X(:)-D(:)).^2))
