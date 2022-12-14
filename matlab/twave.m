% version Mainz
%
%Usage: 
%[h0, h1, AN, BN] = twave (mu, TD, freq, mu_ambient, nth_con, udf_args);
%
% h0 = extinct rate from surface Temperature to back side temperature
% h1 = extinct rate from surface flux to back side temperature
%
% nth_con = [0, 0, ..., 2, 0] (example),
%         0 -> contact
%         1 -> contact with resistance
%              <external boundary>
%                convectCoef = udf_args{pf}(1)
%              <internal boundary>
%                resistance = udf_args{pf}(1)
%         2 -> radiation heat exchange
%              <external boundary>
%                epsR = udf_args{nlayers}(1)
%                tmp0 = udf_args{nlayers}(2)
%              <internal boundary>
%                epsR0=udf_args{pf}(1)
%                epsR1=udf_args{pf}(2)
%                tmp0=udf_args{pf}(3)
%                tmp1=udf_args{pf}(4)
%                (optional) L0=udf_args{pf}(5)
%                (optional) L1=udf_args{pf}(6)
%
%
%Effusity = sqrt(C*rho*lambda), which is intrinsic property of the materials, e.g. Effusity = 22455 for Aluminum plate, more detail please refers https://thermaleffusivity.com/thermaleffusivityvsthermaleffusance/
%
%ThermalDepth = Thickness / sqrt(kappa), where kappa is the diffusion coeffiecent, e.g. kappa = 8.17e-5 for Aluminum plate.
%
%Frequency is the frequency of the input thermal signal, e.g. Frequency = 1e-4 for 0.1mHz disturbance.
%
%Effusity_ambient refers to the effusity of the backside media, e.g. Effusity_ambient = 0 for adiabatic condition.
%
%mu_ambient=4*(5.67e-8)*(epsilon_d)*(T_d)^3/sqrt(2*pi*freq)*exp(-pi/4*i), for vacuum radiation condition.
%
%Extinct as output represents the complex amplification attenuation from the input signal at frontside to the temperature fluctuation at backside.

function [h0,h1,AN,BN]=twave(data_mu,data_TD,freq,mu_ambient,nth_con,udf_args)
if nargin<4
mu_ambient=0;
end
if nargin<5
nth_con=zeros(length(data_TD),1);
%udf_args{1:length(data_TD)}=zeros(1,3);
end
if length(data_mu)~=length(data_TD)
disp('incompatibale input');
return;
end
nlayers=length(data_mu);
data_mu=reshape(data_mu,1,nlayers);
data_TD=reshape(data_TD,1,nlayers);

SRF=sqrt(pi*freq);
An=1;
Bn=0;
AN=ones(nlayers+1,1);
BN=zeros(nlayers+1,1);

switch(nth_con(nlayers))
case 0
	fprintf("manual mode mu_a=%f\n",mu_ambient);
case 1
	cvtCoef=udf_args{nlayers}(1);
	mu_ambient=exp(i*pi/4)*cvtCoef./sqrt(2*pi*freq);
case 2
	epsR=udf_args{nlayers}(1);
	tmp0=udf_args{nlayers}(2);
	mu_ambient=exp(i*pi/4)*4*(5.67e-8)*epsR*tmp0^3./sqrt(2*pi*freq);
otherwise
	fprintf("unknown boundary type at nth_con(%d)\n",nlayers);
	return;
end

newAn=exp((1-i)*SRF*data_TD(nlayers))/2.*((1+mu_ambient/data_mu(nlayers)).*An+(1-mu_ambient/data_mu(nlayers)).*Bn);
newBn=exp(-(1-i)*SRF*data_TD(nlayers))/2.*((1-mu_ambient/data_mu(nlayers)).*An+(1+mu_ambient/data_mu(nlayers)).*Bn);

An=newAn;Bn=newBn;

AN(nlayers)=An;BN(nlayers)=Bn;

for pf=nlayers-1:-1:1
switch(nth_con(pf))
case 0
	newAn=exp((1-i)*SRF*data_TD(pf))/2.*((1+data_mu(pf+1)/data_mu(pf))*An+(1-data_mu(pf+1)/data_mu(pf))*Bn);
	newBn=exp(-(1-i)*SRF*data_TD(pf))/2.*((1-data_mu(pf+1)/data_mu(pf))*An+(1+data_mu(pf+1)/data_mu(pf))*Bn);

case 1
	Rc=udf_args{pf}(1);
	qn=data_mu(pf+1)*exp(-i*pi/4)*sqrt(2*pi*freq).*(An-Bn);
	newAn=exp((1-i)*SRF*data_TD(pf))/2.*(Rc*qn+(1+data_mu(pf+1)/data_mu(pf))*An+(1-data_mu(pf+1)/data_mu(pf))*Bn);
	newBn=exp(-(1-i)*SRF*data_TD(pf))/2.*(Rc*qn+(1-data_mu(pf+1)/data_mu(pf))*An+(1+data_mu(pf+1)/data_mu(pf))*Bn);
case 2
	epsR0=udf_args{pf}(1);
	epsR1=udf_args{pf}(2);
	T0=udf_args{pf}(3);
	T1=udf_args{pf}(4);
	mu0=data_mu(pf);
	mu1=data_mu(pf+1);
	if(length(udf_args{pf})>4)
		L0=udf_args{pf}(5);
		L1=udf_args{pf}(6);
	else
		L0=0;L1=0;
	end
	newAn=exp((1-i)*SRF*data_TD(pf))/2.*((mu1*sqrt(2*pi*freq)*exp(-pi/4*i)/(4*(5.67e-8)*epsR1*T0^3) + (epsR0+L0)*mu1/(epsR1*mu0)) * (An-Bn) + ((1+L1/epsR1)*(T1/T0)^3 + 4*(5.67e-8)*(L0+L1*epsR0/epsR1+L0*L1/epsR1)*T1^3 / (mu0*sqrt(2*pi*freq)*exp(-pi/4*i))) * (An+Bn));
	newBn=exp(-(1-i)*SRF*data_TD(pf))/2.*((mu1*sqrt(2*pi*freq)*exp(-pi/4*i)/(4*(5.67e-8)*epsR1*T0^3) - (epsR0+L0)*mu1/epsR1/mu0) * (An-Bn) + ((1+L1/epsR1)*(T1/T0)^3 - 4*(5.67e-8)*(L0+L1*epsR0/epsR1+L0*L1/epsR1)*T1^3 / (mu0*sqrt(2*pi*freq)*exp(-pi/4*i))) * (An+Bn));

otherwise
	fprintf("unknown boundary type at nth_con(%d)\n",pf);
	return;
end
An=newAn;Bn=newBn;
AN(pf)=An;BN(pf)=Bn;
end

h0=1./(An+Bn);
h1=1./(data_mu(1)*exp(-pi/4*i)*sqrt(2*pi*freq).*(An-Bn));
end
