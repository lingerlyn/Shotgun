%% get rates from calcium - works!!
delta_function_array=zeros(T,1); %does not have to be T  long - just longer then the filter
delta_function_array(end/2)=1;

Gx = @(x,mode) G_mat(x,mode,T,params.calcium_model.g,0);
h_impulse=Gx(delta_function_array,1);
h_tag_impulse=Gx(delta_function_array,1);
h=sum(h_impulse,1);

m=(mean(Y,2)-params.calcium_model.baseline)/h;

%% get Cxx and Cxy from calcium
% N=6;
% Pyy=[]; Pss=[];
% 
% for kk=1:N
%     [Pyy(:,end+1),w]=periodogram(Y(kk,:)+eps*i);
%     Pss(:,end+1)=(Pyy(:,kk)-params.calcium_model.sn^2)./abs(H_impulse.*H_tag_impulse).^2;
% end
% 
% H_impulse=fft(h_impulse,length(w));
% H_tag_impulse=fft(h_tag_impulse,length(w));
% % semilogy(w,Pyy(:,1)-params.calcium_model.sn^2,w,1./abs(H_impulse.*H_tag_impulse).^2);
% semilogy(w,Pyy(:,1))

%% Get S from Y
g=params.calcium_model.g;
sn=params.calcium_model.sn;
h_inv_impulse=filter([1;-g(:)],1,delta_function_array);
h_inv2_impulse=filter([1;-g(:)],1,flipud(h_inv_impulse));
h_inv2_sum=sum(h_inv2_impulse);

% deconvlution
% zeros=[1;-g(:)];poles=1; % inverse filter - bad perfomance on covariance
for kk=1:N
    Sxx=m(kk)*(1-m(kk));
    sn_normalized=sn/sqrt(Sxx);
    a=sn_normalized^2;
    ze=[0;-1;g]; po=[g*a; -(g^2*a+a+1) ; g*a]; % MinMax filter - unstable
    spikes(kk,:)= filter(ze,poles,Y(kk,:)-params.calcium_model.baseline);
end
mean(spikes,2)
mean(true_spikes,2)
