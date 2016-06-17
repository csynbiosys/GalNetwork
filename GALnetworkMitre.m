% [T,Y]=ode45(@GALnetworkMitre,[0 1000],ones(5,1))
function dy = GALnetworkMitre(t,y)
t
%Definition of the parameters
ktr3 = 0.329;
ktr80 = 0.147;
ktr2 = 0.678;
ktr1 = 1.042;
ktl3 = 645;
ktl80 = 210;
ktl2 = 800;
ktl1 = 187;
% c = 4.215*e+7;
ua = 4.438e-3;
ym = 4.332e-2; %There should be a value for each gene but I can find the other ones
yg3 = 7.112e-3;
yg80 = 2.493e-3;
% yg2 = 0;
yg1 = 0;
kd80 = 3e-7;
kb80 = 5e-6;
kb3 = 6e-8;
kb1 = 6e-8;
kd3 = 1.25e-2;
kd1 = 1;
ks = 4000;
kc3 = 0.5;
kc1 = 8e-5;
a = 4350;
kgk = 702;
d = 59200;
K = 1;
km = 1.2;
kic = 160;
kiu = 19.1;
ub = 0.00512;
uc = 0.3611;
nu = 1;
yb = 0.001416;
yc = 0.8592;
ny = 1;
xc = 0.2443;
nx = 1;
yb = 0.0003;
yc = 2.9989;
k80 = sqrt(kd80*kb80);
k3 = ((sqrt(kd3*kb3)*(yg3+ua))/kc3);
k1 = ((sqrt(kd1*kb3*kb1)*(yg1+ua))/kc1);
Ge =  0.03;  % External galactose in w/v %
R  =  0;  % Glucose 
syms z; % creating a symbolic variable z

R1 = 1-(1/(1+((k80/y(2))^2)+((y(1)*y(5)/(k3*ks+k3*y(5)))^2)+((y(4)*y(5)/(k1*ks+k1*y(5)))^2)));
R4 = 1-(1/(1+symsum(((k80/y(2))^(2*z)),z,1,4)+symsum(((y(1)*y(5)/(k3*ks+k3*y(5)))^(2*z)),z,1,4)+symsum(((y(4)*y(5)/(k1*ks+k1*y(5)))^(2*z)),z,1,4)));
u = ua+((ub*R^nu)/((uc^nu)+(R^nu)));
yg2 = (yb*(R^ny))/((yc^ny)+(R^ny));
x = 1/(((R/xc)^nx)+1);
yr = (1-yb)+(yb/(yc+R));
dy = zeros(5,1);
% Phosphorylation function
P = (kgk*kiu*kic/(kiu*km+kic*y(5)));
Kp = (km+y(5))*kiu*kic/((kiu*km)+(kic*y(5)));

%% Differential equations 

% Gal3p
dy(1) = (ktl3*ktr3*x/(ym+u))*R1-y(1)*(yg3+u+(kc3*y(5)/(ks+y(5))));
% Gal80p
dy(2) = (ktl80*ktr80/(ym+u))*R1-y(2)*(yg80+u);
% Gal2p
dy(3) = (ktl2*ktr2/(ym+u))*x*R4-y(3)*(yg2+u);
% Gal1p
dy(4) = (ktl1*ktr1/(ym+u))*x-y(4)*(yg1+u+(kc1*y(5)/(ks+y(5))));
% Internal galactose 
dy(5) = a*yr*y(3)*((Ge/(K+Ge))-(y(5)/(K+y(5))))-(2*P*y(5)*y(4)/(Kp+sqrt((Kp^2)+(4*P*y(5)*y(4)/d))))-y(5)*((kc3*y(5)/(ks+y(5)))+(kc1*y(5)/(ks+y(5))))-u*y(5);


end 