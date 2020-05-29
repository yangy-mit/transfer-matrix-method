% Function handles for Transfer Matrix Method (TMM) to perform
% electromagnetic analysis of multilayer media
% 
% Both TE and TM polarizations are considered
% 
% The refractive index can be complex (materials with loss or gain)
% 
% Developed by: Yujia Yang, MIT, 2013
% Updated by: Yujia Yang, MIT, 2020

%Constants
c_0=299792458; 
eps_0=8.85e-12;
mu_0=4*pi*1e-7;

%Wavevector in propagation direction, kz
kz = @(omega,n,angle) omega.*n./c_0.*cos(angle);

%Impedance Z0
Z0 = @(n,angle,pol) pol.*sqrt(mu_0./eps_0/n.^2)./cos(angle)+(1-pol).*sqrt(mu_0./eps_0/n.^2).*cos(angle);
%TE:pol=1;TM:pol=0

%p10
p10 = @(n1,n0,angle1,angle0,pol) pol.*Z0(n1,angle1,pol)./Z0(n0,angle0,pol)+(1-pol).*Z0(n0,angle0,pol)./Z0(n1,angle1,pol);

%Reflection Coefficient gamma10
gamma10 = @(n1,n0,angle1,angle0,pol) (1-p10(n1,n0,angle1,angle0,pol))./(1+p10(n1,n0,angle1,angle0,pol));

%Transfer Matrix
TM = @(omega,n1,n0,thick1,angle1,angle0,pol) 0.5.*(1+p10(n1,n0,angle1,angle0,pol))*[exp(1i*kz(omega,n1,angle1).*thick1),gamma10(n1,n0,angle1,angle0,pol).*exp(1i*kz(omega,n1,angle1).*thick1); ...
                                                                                    gamma10(n1,n0,angle1,angle0,pol).*exp(-1i*kz(omega,n1,angle1).*thick1),exp(-1i*kz(omega,n1,angle1).*thick1)];
