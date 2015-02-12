function [g,sn,z,b]= armafit(Y,order,varargin)
% input:
% Y - calcium trace
% order - AR model order
% output:
% g - AR model rates
% sn - noise std (in relevant format for inference code)
% z - MA model rates
% b - bias
% varargin - previously identified spikes and bias (optional)
 
if isempty(varargin)  % sysid without input
    dat = iddata(Y,Y*0+1);
    sys=polyest(dat,[0,1,order,order,0,0]);
    g=-sys.D(2:end)';
    sz=sys.C(end)/(sys.D(end)-sys.C(end));
    z=[1 (sys.C(2:end)*(sz+1)-sz*sys.D(2:end))];
    b=sys.B(1);
    sn=(sz/(sz+1))*sqrt(sys.NoiseVariance/2);
else %iterated sysid
    U=varargin{1}; %spikes and bias
    spikes=U(1:end-1)';
    bias=U(end);
    dat = iddata(Y-bias,spikes);
%     sys=polyest(dat,[order,order,order,0,0,0]);
%     g=-sys.A(2:end)';    
    sys=polyest(dat,[0,order,order,order,order,0]);
    g=-sys.F(2:end)';    
    z=sys.B;    
    sn=sqrt(sys.NoiseVariance/2);
    b=bias;
end
    
end

