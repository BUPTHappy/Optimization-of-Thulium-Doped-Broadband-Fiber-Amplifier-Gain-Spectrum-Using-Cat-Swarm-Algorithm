
clear;
GL=zeros(8,1);
% first curve
L=2; % fiber length (m)
Ps_0=1.0e-6;  % signal power at input end: -30dBm
Pp_0=200e-3; % 980nm main pump power at input end: 200mW
Pase_0=0;     % ASE in signal bandwidth power at input end: 0
N00=1*1.6e25;

x=1450:1:1520;
y=normpdf(x,1460,255);
temp1=(y/max(y))*(1.0e-24);  %取y最大值



sig_se=temp1;
sig_sa=temp1;
swl=x*1e-9;
c=3e8;
h=6.626e-34; % the Plank constant

len=length(swl);
G=zeros(1,len);
for count=1:len
    [z,P]=ode45(@fun,[0,L],[Ps_0,Pp_0,Pase_0],'AbsTol',N00,swl(count),sig_se(count),sig_sa(count));
    G1(1,count)=10*log10(P(length(P),1)./Ps_0);
    Ps(1,count)=P(length(P),1);
    Pase(1,count)=P(length(P),3);
end


N11=0.4*1e25;% fiber length (m)
len=length(swl);
G2=zeros(1,len);
for count=1:len
    [z,P]=ode45(@fun,[0,L],[Ps_0;Pp_0;Pase_0],'AbsTol',N11,swl(count),sig_se(count),sig_sa(count));
    G2(1,count)=10*log10(P(length(P),1)./Ps_0);
    Ps(1,count)=P(length(P),1);
    Pase(1,count)=P(length(P),3);
end


N22=0.6*1e25; % fiber length (m)
len=length(swl);
G3=zeros(1,len);
for count=1:len
    [z,P]=ode45(@fun,[0,L],[Ps_0;Pp_0;Pase_0],'AbsTol',N22,swl(count),sig_se(count),sig_sa(count));
    G3(1,count)=10*log10(P(length(P),1)./Ps_0);
    G33(1,count)=P(length(P),1)./Ps_0;
    Ps(1,count)=P(length(P),1);
    Pase(1,count)=P(length(P),3);
end


N33=0.8*1e25; % fiber length (m)
len=length(swl);
G4=zeros(1,len);
for count=1:len
    [z,P]=ode45(@fun,[0,L],[Ps_0;Pp_0;Pase_0],'AbsTol',N33,swl(count),sig_se(count),sig_sa(count));
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



