function dy = montiGALphilipp(t,y)
% [T,Y]=ode45(@montiGALmodel,[0 5000],ones(4,1))
% This gal network model is from the paper "A modified galactose network
% with implication for growth"
%t
%Defining the model parameters
kf81 = 100;
kr81 = 1500;
kf83 = 100;
kr83 = 1;
kf84 = 100;
kr84 = 25;
aG1 = 15;
a0G1 = 0.418;
aG3 = 0.9;
aG4 = 0.2;
a0G80 = 0.6;
aG80 = 0.9;
KG1 = 8;
KG3 = 8;
KG80 = 2;
n1 = 3;
n3 = 2;
n80 = 2;
yG1 = 0.004;
yG3 = 0.004;
yG4 = 0.004;
yG80 = 0.004;
yC81 = 0.004;
yC83 = 0.004;
yC84 = 0.004;
e = 0.1;
C = 22.1;
KM = 0.086;
% s=10E-3*(sign(t-1000)-sign(2000-t)+1)/2; %Input step after 1000 for 1000 minutes 
s=10E-4*(sign(t-1000)+1)/2-10E-4*(sign(t-2000)+1)/2;



as = C*s/(s+KM);
w = (kr81*kf81/(kr81+yC81))-kf81;
d = (kr83*kf83/(kr83+yC83))-kf83;
b = (kr84*kf84/(kr84+yC84))-kf84;

dy = zeros(4,1);
dy(1) = a0G1 + as*e+aG1*((y(3)^n1)/((y(3)^n1)+ KG1^n1))+ w*y(1)*y(4)-yG1*y(1);
dy(2) = as + aG3*((y(3)^n3)/(y(3)^n3+KG3^n3))+d*y(2)*y(4)-yG3*y(2);
dy(3) = aG4 + b*y(3)*y(4)-yG4*y(3);
dy(4) = a0G80 + aG80*((y(3)^n80)/((y(3)^n80)+(KG80^n80)))+w*y(1)*y(4)+b*y(3)*y(4)+d*y(2)*y(4)-yG80*y(4);



end


