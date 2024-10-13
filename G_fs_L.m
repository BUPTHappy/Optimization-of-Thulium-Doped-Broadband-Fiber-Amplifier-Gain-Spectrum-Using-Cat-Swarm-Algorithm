
clear;
GL=zeros(8,1);
% first curve
L1=2; % fiber length (m)
Ps_0=1.0e-6;  % signal power at input end: -30dBm
Pp_0=200e-3; % 980nm main pump power at input end: 200mW
Pase_0=0;     % ASE in signal bandwidth power at input end: 0
N=1*1.6e25;


x=1450:1:1520;
a0=0.2483; 
a1=0.1299;
b1=0.0409;
a2=0.02744;
b2=0.1377;
a3=0.05425;
b3=0.0486;
a4=-0.0163;
b4=-0.06018 ;
a5=-0.007029;
b5=0.008691;
a6=0.03024;
b6=0.0005771;
w=0.007206;
y1=a0+a1*cos(x*w)+b1*sin(x*w)+a2*cos(2*x*w)+b2*sin(2*x*w)+a3*cos(3*x*w)+b3*sin(3*x*w)+a4*cos(4*x*w)+b4*sin(4*x*w)+a5*cos(5*x*w)+b5*sin(5*x*w)+a6*cos(6*x*w)+b6*sin(6*x*w);
temp1=(y1/max(y1))*(1.0e-24);  %取y最大值

%吸收截面与波长关系拟合
A0=0.6494;
A1 =0.1692;
B1=0.0605;
A2=0.05075;
B2=0.05163;
A3=0.004536;
B3=0.01396;
w=0.05249;
y2=A0+A1*cos(x*w)+B1*sin(x*w)+A2*cos(2*x*w)+B2*sin(2*x*w)+A3*cos(3*x*w)+B3*sin(3*x*w);
temp2=(y2/max(y2))*(1.0e-24);  %取y最大值


sig_se=temp1;
sig_sa=temp2;
swl=x;


len=length(swl);
G=zeros(1,len);
for count=1:len
    [z,P]=ode45(@fun,[0,L1],[Ps_0;Pp_0;Pase_0],'AbsTol',N,swl(count),sig_se(count),sig_sa(count));
    G1(1,count)=10*log10(P(length(P),1)./Ps_0);
    Ps(1,count)=P(length(P),1);
    Pase(1,count)=P(length(P),3);
end



L2=4; % fiber length (m)
len=length(swl);
G=zeros(1,len);
for count=1:len
    [z,P]=ode45(@fun,[0,L2],[Ps_0;Pp_0;Pase_0],'AbsTol',N,swl(count),sig_se(count),sig_sa(count));
    G2(1,count)=10*log10(P(length(P),1)./Ps_0);
    Ps(1,count)=P(length(P),1);
    Pase(1,count)=P(length(P),3);
end



L3=6; % fiber length (m)
len=length(swl);
G=zeros(1,len);
for count=1:len
    [z,P]=ode45(@fun,[0,L3],[Ps_0;Pp_0;Pase_0],'AbsTol',N,swl(count),sig_se(count),sig_sa(count));
    G3(1,count)=10*log10(P(length(P),1)./Ps_0);
    Ps(1,count)=P(length(P),1);
    Pase(1,count)=P(length(P),3);
end


L4=8; % fiber length (m)
len=length(swl);
G=zeros(1,len);
for count=1:len
    [z,P]=ode45(@fun,[0,L4],[Ps_0;Pp_0;Pase_0],'AbsTol',N,swl(count),sig_se(count),sig_sa(count));
    G4(1,count)=10*log10(P(length(P),1)./Ps_0);
    Ps(1,count)=P(length(P),1);
    Pase(1,count)=P(length(P),3);
end

figure(1);
plot(swl,G1,'r.');
text(1.61e-6,16,' \leftarrow L=2','FontSize',11)
grid on;
xlabel('Wavelength (nm)');
ylabel('Gain (dB)');
hold on;

plot(swl,G2,'b.-');
plot(swl,G3,'g');
plot(swl,G4,'m--');



