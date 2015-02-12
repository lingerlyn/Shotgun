function P = GetParams(Y,order,noise_method,timescale_method,varargin)
% estimate a global time constants, bias and noise power for each pixel
% input:
% Y - NxT calcium trace
% order - order of AR filter
% noise_method - method used to infer noise power. Can be 'keep','psd' or 'mcmc'
% timescale_method - method used to infer input timescales. Can be 'arpfit' or 'sysid'
% varargin - Nx(T+1) array of spikes and bias (in t=T+1)
% output:
% P - parameter struct of size N:
% P.g - poles polynomial coeffieints 
% P.Cb - bias
% P.sn - noise stdev
% P.z - zeros

[N,T] = size(Y);
% P=cell(N,1);

for nn=1:N
    switch timescale_method        
        case 'arpfit'
            P{nn} = arpfit(Y(nn,:),order);
            P{nn}.z=[];
        case 'sysid'
            if isempty(varargin)
                [g,sn,z,Cb]= armafit(Y(nn,:)',order);
            else %iterated sysid
                input=varargin{1};                   
                [g,sn,z,Cb]= armafit(Y(nn,:)',order,input(nn,:));
            end
            P{nn}.g=g;
            P{nn}.sn=sn;
            P{nn}.z=z;
            P{nn}.Cb=Cb;
        otherwise
            error('unknown timescale method')
    end
    
    switch noise_method
        case 'keep'
            continue
        case 'psd'
            [ sn, ~ ] = GetSnPSD(Y(nn,:)');
            P{nn}.sn=sn;
        case 'mcmc'
            error('mcmc noise method not written yet - add this!')
        otherwise
            error('unknown noise method')
    end
end

end

