%fun.m
% the basic function which is called by the function ode45()
function f=fun(z,P,N,swl,sig_se,sig_sa)
% P=[Ps,Pp,Pase];
% NB=N;
c=3e8; % the velocity of lightwave (m/s)
h=6.626e-34; % the Plank constant
alpha=0.1; % background loss (dB/m)

%---------------------------------------------------------------
% Spontaneous emission rate (1/s)

A32=1;
A21=108.6;
%---------------------------------------------------------------
% Cross-sections (m2)
sigma_se=sig_se; % stimulated emission cross-section of signal
sigma_sa=sig_sa; % absorption cross-section of signal
sigma_23ase_a=sigma_sa;
sigma_23ase_e=sigma_se;
sigma_13pe=0.8e-24;  %stimulated emission cross-section of pump
sigma_13pa=2*3.5e-24; %absorption cross-section of pump

%---------------------------------------------------------------


%---------------------------------------------------------------
r=3.0e-6; % core radius of the fiber (m)
a=r;
lambda_P=980e-9; % wavelength of main pump
lambda_s=swl*1e-9;%wavelength of signal
lambda_ase=lambda_s;% wavelength of  ASE

f_p=c/lambda_P; % frequency of main pump
f_s=c/lambda_s; % frequency of signal
f_ase=c/lambda_ase; % frequency of 800nm ASE

delta=c*1e-9/lambda_s^2;% ASE frequency bandwidth

% Overlapping factors
v_p=0.22*2*pi*a/lambda_P;
r_p=a/sqrt(2)*(0.65+1.619*v_p^(-1.5)+2.879*v_p^(-6));
Tp=1-exp(-a^2/r_p^2);

v_s=0.22*2*pi*a/lambda_s;
r_s=a/sqrt(2)*(0.65+1.619*v_s^(-1.5)+2.879*v_s^(-6));
Ts=1-exp(-a^2/r_s^2);
Tase=Ts;

%---------------------------------------------------------------
% P=[Ps;Pp;Pase] is a column vector
Aeff=pi*r^2;

W12sa=sigma_sa*P(1)/(h*f_s*Aeff); % W12sa
W12se=sigma_se*P(1)/(h*f_s*Aeff); % W12se
Wp=sigma_13pa*P(2)/(h*f_p*Aeff); % W12pa
W13pe=sigma_13pe*P(2)/(h*f_p*Aeff); % W12pe
W12ase_a=sigma_23ase_a*P(3)/(h*f_ase*Aeff); % W12ase_a
W12ase_e=sigma_23ase_e*P(3)/(h*f_ase*Aeff); % W12ase_e

W23=W12sa+W12ase_a;
W32=W12se+W12ase_e;
%---------------------------------------------------------------
% N0, N1, N3
N1=(A21*A32*N + A21*N*W32)/(W23*Wp + W32*Wp + A21*A32 + A21*W32 + A21*Wp + A32*Wp);
N2=(N*Wp*(A32 + W32))/(W23*Wp + W32*Wp + A21*A32 + A21*W32 + A21*Wp + A32*Wp);
N3=(N*Wp*(A21 + W23))/(W23*Wp + W32*Wp + A21*A32 + A21*W32 + A21*Wp + A32*Wp);

%---------------------------------------------------------------
% Note:the column vector P=[Ps;Pp;Pase]
f=[ P(1)*Ts*(sigma_se*N3-sigma_sa*N2)-alpha*P(1);
    P(2)*Tp*(sigma_13pe*N3-sigma_13pa*N1)-alpha*P(2);   
    P(3)*Tase*(sigma_se*N3-sigma_sa*N2)+2*sigma_se*N3*Ts*h*f_s*delta-alpha*P(3)];
