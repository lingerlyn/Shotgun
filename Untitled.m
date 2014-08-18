% a= COV(1:N,N+1:end);%%%%%;
% b=COV_res(1:N,N+1:end);%%%%%;
% plot(a(:),b(:))
% 
% a=COV_res(1:N,1:N)
% b=COV(1:N,1:N)
% plot(a(:),b(:))

COV=[Cxx Cxy; Cxy' Cxx];
z=eig(COV);
sum(z<0)
% scatter(real(z),imag(z));