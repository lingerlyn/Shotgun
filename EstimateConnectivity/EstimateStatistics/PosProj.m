function [CXX,CXY ] = PosProj( CXX,CXY )
% Do a positive semidefinite projection

COV = [CXX CXY; CXY' CXX];
N=size(CXX,1);

if(any(eig(COV)<0))
    disp('COV is not positive semidefinite;')
        disp('correcting...');
        [v,d]=eig(COV);
        COV=v*spdiags(max(diag(d),0),0,2*N,2*N)*v';            
end

CXY = COV(1:N,N+1:end);
CXX = COV(1:N,1:N);

end

