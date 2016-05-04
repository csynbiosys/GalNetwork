function dy = montiGALmodel(t,y)
s=10E-4*(sign(t-1000)+1)/2-10E-4*(sign(t-2000)+1)/2;
% Parameters
kf81=100;
kr81=1500;
kf83=100;
kr83=1;
kf84=100;
kr84=25;
ag1=15;
a0g1=0.418;
ag3=0.9;
ag4=0.2;
a0g80=0.6;
%apg80=0.994;
ag80=0.9;
kg1=8;
kg3=8;
kg80=2;
n1=3;
n3=2;
n80=2;
gammag1=0.004;
gammag3=0.004;
gammag4=0.004;
gammag80=0.004;
gammac81=0.004;
gammac83=0.004;
gammac84=0.004;
epsilon=0.1;
C=22.1;

%Those values are actually negative
w = (kr81*kf81/(kr81+yC81))-kf81;
d = (kr83*kf83/(kr83+yC83))-kf83;
b = (kr84*kf84/(kr84+yC84))-kf84;

dy = zeros(4,1);
dy(1) = a0G1 + as*e+aG1*((y(3)^n1)/((y(3)^n1)+ KG1^n1))+ w*y(1)*y(4)-yG1*y(1);
dy(2) = as + aG3*((y(3)^n3)/(y(3)^n3+KG3^n3))+d*y(2)*y(4)-yG3*y(2);
dy(3) = aG4 + b*y(3)*y(4)-yG4*y(3);
dy(4) = a0G80 + aG80*((y(3)^n80)/((y(3)^n80)+(KG80^n80)))+w*y(1)*y(4)+b*y(3)*y(4)+d*y(2)*y(4)-yG80*y(4);

end 