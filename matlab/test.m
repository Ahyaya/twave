mu=sqrt([11.5*50*945,5.4e-3*64*1090,11.5*50*945]);
td=[0.010./sqrt(11.5/50/945),0.02./sqrt(5.4e-3/64/1090),0.020./sqrt(11.5/50/945)];
con=[0,2,2];
udf{2}(1)=0.026;
udf{2}(2)=0.026;
udf{2}(3)=321.90;
udf{2}(4)=294.83;
udf{2}(5)=9.57e-3;
udf{2}(6)=9.57e-3;
udf{3}(1)=0.05;
udf{3}(2)=293.15;

freq=1e-4;
Ts=360.62;
[h0,h1,An,Bn]=twave(mu,td,freq,0,con,udf);
hs=0.64/(1/h1+4*(5.67e-8)*0.79*Ts^3/h0);

%format long;
%disp([An,Bn]);

fprintf("h0(f) = %f\n",abs(h0));
fprintf("h1(f) = %f\n",abs(h1));
fprintf("hs(f) = %f\n",abs(hs));
dI=0.175*freq^(-1/3);
fprintf("\nSCF induced fluctuation:\n%f %f K/sqrt{Hz}\n",dI*abs(hs)/abs(h0),dI*abs(hs));
