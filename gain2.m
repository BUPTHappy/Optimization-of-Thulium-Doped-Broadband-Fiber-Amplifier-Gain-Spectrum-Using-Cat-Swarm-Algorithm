function y=gain2(x)
L1=x(1); % fiber length (m) 设为其中一个优化变量
Ps_0=1.0e-6;  % signal power at input end: -30dBm
Pp_0=200e-3; % 980nm main pump power at input end: 200mW
Pase_0=0;     % ASE in signal bandwidth power at input end: 0
N=x(2)*1e25; %设为另一个优化变量

xx=1450:10:1520;
y1=normpdf(xx,1480,60);
y2=normpdf(xx,1500,60);
temp1=(y1/max(y1))*(1.0e-25); 
temp2=(y2/max(y2))*(1.0e-25);

% xx=1450:1:1520;
% a0=0.2483; 
% a1=0.1299;
% b1=0.0409;
% a2=0.02744;
% b2=0.1377;
% a3=0.05425;
% b3=0.0486;
% a4=-0.0163;
% b4=-0.06018 ;
% a5=-0.007029;
% b5=0.008691;
% a6=0.03024;
% b6=0.0005771;
% w=0.007206;
% y1=a0+a1*cos(xx*w)+b1*sin(xx*w)+a2*cos(2*xx*w)+b2*sin(2*xx*w)+a3*cos(3*xx*w)+b3*sin(3*xx*w)+a4*cos(4*xx*w)+b4*sin(4*xx*w)+a5*cos(5*xx*w)+b5*sin(5*xx*w)+a6*cos(6*xx*w)+b6*sin(6*xx*w);
% temp1=(y1/max(y1))*(1.0e-25);  %取y最大值
% 
% %吸收截面与波长关系拟合
% A0=0.6494;
% A1 =0.1692;
% B1=0.0605;
% A2=0.05075;
% B2=0.05163;
% A3=0.004536;
% B3=0.01396;
% w=0.05249;
% y2=A0+A1*cos(xx*w)+B1*sin(xx*w)+A2*cos(2*xx*w)+B2*sin(2*xx*w)+A3*cos(3*xx*w)+B3*sin(3*xx*w);
% temp2=(y2/max(y2))*(1.0e-25);  %取y最大值


sig_se=temp1;
sig_sa=temp2;
swl=xx;


len=length(swl);
G=zeros(1,len);
for count=1:len %将L1设置成第四行设置的光纤长度变量，将N改为第八行设置的掺杂浓度变量
    [z,P]=ode45(@fun,[0,L1],[Ps_0;Pp_0;Pase_0],'AbsTol',N,swl(count),sig_se(count),sig_sa(count));
    G1(1,count)=10*log10(P(length(P),1)./Ps_0);
    Ps(1,count)=P(length(P),1);
    Pase(1,count)=P(length(P),3);
end
y=sum(G1);
end




