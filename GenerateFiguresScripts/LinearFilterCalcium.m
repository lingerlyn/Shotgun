g=params.calcium_model.g;
sn=params.calcium_model.sn;
delta_function_array=zeros(T,1); %does not have to be T  long - just longer then the filter
delta_function_array(1)=1;

Gx = @(x,mode) G_mat(x,mode,T,g,0);
h_impulse=Gx(delta_function_array,1);
h_tag_impulse=Gx(delta_function_array,1);

rates=(mean(Y,2)-params.calcium_model.baseline)/sum(h_impulse,1);

for kk=1:N
% Calculate power spectral density of input image
Sxx=rates(kk)*(1-rates(kk));

% Calculate power spectral density of input image
Snn = sn.^2;

% Calculate Fourier Transform of LSI Filter
H = fft(h_impulse);
H2 = abs(H).^2;

% Calculate thresholded 1/H
HINV = H.^(-1);
index = find(abs(H) < 0.8);
hzeros = length(index);    % Return number of elements below threshold
HINV(index) = 0;

% Calculate Wiener Filter
F=(H2.*Sxx)./((H2.*Sxx) + Snn); %denoising part
F(index) = 0; %just for the correct calculation of the normalization constants later
G = HINV.* F;

% Restore Image
Y_f = fft(Y(kk,:)-params.calcium_model.baseline);
S_f = Y_f.*G';
spikes_est(kk,:) = ifft2(S_f);
end


% 
for kk=1:N
    subplot(N,1,kk)
    hold off
    plot(spikes_est(kk,:))
    hold all
    plot(true_spikes(kk,:),'.')
    set(gca,'xtick',[])
    xlim([100 250])
end
    hold off
    
%%
figure(1)
clf(gcf)
imagesc(spikes_est)
xlim([100 250])
figure(2)
clf(gcf)
imagesc(true_spikes)
xlim([100 250])

%%
% rates=mean(spikes_est,2);
rates=(mean(Y,2)-params.calcium_model.baseline)/sum(h_impulse,1);
Cxx=spikes_est*spikes_est'/T-rates*rates';
% Cxx(eye(N)>0.5)=Cxx(eye(N)>0.5)-sn^2*abs(mean(abs(G).^2));
Cxx(eye(N)>0.5)=rates.*(1-rates);
Cxy=spikes_est(:,1:end-1)*spikes_est(:,2:end)'/T-rates*rates';
w=(1:T)*2*pi/T;
Cxy(eye(N)>0.5)=Cxy(eye(N)>0.5)+sn^2*abs(mean(exp(1i*w').*abs(G).^2));
Cxy=Cxy/(mean(abs(F).^2));
%%
rates_true=mean(true_spikes,2); %% get rates from calcium - works!!
% rates_true=mean(true_spikes,2);
Cxx_true=true_spikes*true_spikes'/T-rates_true*rates_true';
Cxy_true=true_spikes(:,1:end-1)*true_spikes(:,2:end)'/T-rates_true*rates_true';

%%
mi=min(Cxx_true(:));ma=max(Cxx_true(:));
figure(1)
clf(gcf)
imagesc(Cxx,[mi ma])
figure(2)
clf(gcf)
imagesc(Cxx_true,[mi ma])
figure(3)
clf(gcf)
mi=min([Cxx_true(:)])*1.2;
ma=max([Cxx_true(:)])*1.2;
A_ind=linspace(mi,ma,100);
plot(A_ind,A_ind,'k--','linewidth',1);
hold all
plot(Cxx_true(:),Cxx(:),'.')

%%
mi=min(Cxy_true(:));ma=max(Cxy_true(:));
figure(1)
clf(gcf)
colorbar
imagesc(Cxy,[mi ma])
figure(2)
clf(gcf)
colorbar
imagesc(Cxy_true,[mi ma])
figure(3)
clf(gcf)
mi=min([Cxy_true(:)])*1.2;
ma=max([Cxy_true(:)])*1.2;
A_ind=linspace(mi,ma,100);
plot(A_ind,A_ind,'k--','linewidth',1);
hold all
plot(Cxy_true(:),Cxy(:),'.')

