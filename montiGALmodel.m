% [t,y]=ode45(@montiGALmodel,[0 3000],[31.9700    0.0001    0.3630   34.1947]);plot(t,y)
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
apg80=0.994;
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
km=0.086;

omega=((kr81*kf81)/(kr81+gammac81))-kf81;
alphas=C*((s)/(s+km));
delta=((kr83*kf83)/(kr83+gammac83))-kf83;
beta=((kr84*kf84)/(kr84+gammac84))-kf84;

dy = zeros(3,1);
dy(1)= a0g1+alphas*epsilon+ag1*((y(3)^n1)/(y(3)^n1+kg1^n1))+omega*y(1)*y(4)-gammag1*y(1);
dy(2)= alphas + ag3*((y(3)^n3)/(y(3)^n3+kg3^n3))+delta*y(2)*y(4)-gammag3*y(2);
dy(3)= ag4+beta*y(3)*y(4)-gammag4*y(3);
dy(4)=a0g80+ag80*((y(3)^n80)/(y(3)^n80+kg80^n80))+omega*y(1)*y(4)+beta*y(3)*y(4)+delta*y(2)*y(4)-gammag80*y(4);

end 