
function dy = VenturelliModel(t,y)
% [T,Y]=ode45(@VenturelliModel,[0 3000],ones(4,1))
%Parameters Venturelli model
t
kf81 = 100; 
kr81 = 1500;
kf83 = 100;
kr83 = 1;
kf84 = 100;
kr84 = 25;
ag1 = 15;
ag3 = 0.9;
ag4 = 0.2;
a0g80 = 0.6;
ag80 = 0.9;
kg1 = 8;
kg3 = 8;
kg80 = 2;
n1 = 3;
n3 = 2;
n80 = 2;
yg1 = 0.004;
yg3 = 0.004;
yg4 = 0.004;
yg80 = 0.004;
yc81 = 0.004;
yc83 = 0.004;
yc84 = 0.004;
e = 0.1;
ag1s = 0.1;
ag3s = 0.1;
ag80s = 1.5;
agal = 0.4; %galactose at constant input rate in nM/min and values from 0-2

w = (kr81*kf81/(kr81+yc81))-kf81;
d = (kr83*kf83/(kr83+yc83))-kf83;
b = (kr84*kf84/(kr84+yc84))-kf84;


dy = zeros(4,1);
dy(1) = agal*e+((ag1*(y(3))^3)/(((y(3))^3)+kg1^3))+w*y(1)*y(4)-yg1*y(1);
dy(2) = agal+(ag3*(y(3)^2)/(y(3)^2+kg3^2))+d*y(2)*y(4)-yg3*y(2);
dy(3) = ag4+b*y(4)*y(3)-yg4*y(3);
dy(4) = a0g80+(a0g80*(y(3)^2)/((y(3)^2)+kg80^2))+w*y(1)*y(4)+d*y(2)*y(4)+b*y(4)*y(3)-yg80*y(4);
end