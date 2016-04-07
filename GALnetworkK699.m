% [t,y]=ode45(@GALnetworkK699,[0 1000],ones(23,1))
function dy = GALnetworkK699(t,y)
t
% values for galactose network module for the K699-strain
% b1 = 9.92;
% b2 = 6.94;
% b3 = 18.0;
% b4 = 0.86;
% b80= 4.00;
% a1 = 1.09;
% a2 = 1.20;
% a3 = 36.0;
% a80= 3.00;
% d1 = 0.0033;
% y1 = 0.036;
% y2 = 0.026;
% y3 = 0.036;
% y80 = 0.036;
% k3 = 5.0*10^(-8);
% kn3 = 890;
% k4d = 0.10;
% kn4d = 1.0;
% k80d = 0.10;
% kn80d = 170;
% k80 = 0.10;
% kn80 = 0.03;
% kr = 0.10;
% knr = 1.80;
% kcat = 3350;
% kmgk = 1.29*10^7;
% ktr = 4350;
% kmtr = 2.15*10^8;
% atr = 10.0;
% kp = 0.091;
% kq = 0.0556;
% cp = 1;
% cq = 30;
% GALe =2.366*10^8;


% new values from the simplified model
s = 2.366*10^8; %Galactose concentration outside the cell
a0g1 = 106.11;
as = 253.86*22.1*(s/(s+0.086));
e =0.1;
ag1 = 15*253.86;
Kg1 = 8*253.86;
Kg3 = 8*253.86;
Kg80= 2*253.86;
yg1 = 0.004;
yg3 = 0.004;
yg4 = 0.004;
yg80 = 0.004;
a0g80 = 0.6*253.86;
ag80 = 0.9*253.86;
w = ((1500*100*253.860/(1500+0.004))-100*253.86);
d = ((1*100*253.860/(1+0.004))-100*253.86);
b = ((25*100*253.860/(25+0.004))-100*253.86);
ag3 = 0.9*253.86;
ag4 = 0.2*253.86;
n1 = 3;
n3 = 2;
n80 =2; 


% values for glucose network module for the K699-strain

aglu = 215;
eglu = 6.452*10^5;
tglu = 0.4;
dglu = 0.1;
yglu = 0.0633;
dd = 0.0033;
cglu = 1.075*10^7;
b = 1.8;
ktr2 = 4350;
kmtr2 = 6.022*10^8;
atr2 = 1.0;
uglu = 5350;
kglu = 1.29*10^7;
GLUe = 0; %2.957*10^8;
p = 1.29*10^7;
q = 0.8;
% vsd = 9.30*10^(-6);
% ksd = 30;



dy = zeros(7,1);

dy(5)= ktr2*y(7)*((GLUe-y(5))/(kmtr2+GLUe+y(5)+(atr2*GLUe*y(5))/kmtr2))-((uglu*y(5)*y(7))/(kglu+y(5)))-dd*y(5);
dy(5)= (aglu+eglu*((y(5)/cglu)^b))/(1+((y(5)/cglu)^b))-yglu*y(6);
dy(7)= tglu*y(6)-dglu*y(7);
y19= (p^q)/((p^q)+(y(7)^q));
% y20= ktr*y(2)*((GALe-y(15))/(kmtr+GALe+y(15)+(atr*GALe*y(15)/kmtr)));

dy(1) = a0g1+as*e+ag1*((y(3)^n1)/((y(3)^n1)+Kg1^n1))+w*y(1)*y(4)-yg1*y(1);
dy(2) = as + ag3*((y(3)^n3)/((y(3)^n3)+Kg3^n3))+d*y(2)*y(4)-yg3*y(2);
dy(3) = ag4*y19+b*y(3)*y(4)-yg4*y(3);
dy(4) = a0g80+ag80*((y(3)^n80)/((y(3)^n80)+Kg80^n80))+w*y(1)*y(4)+b*y(3)*y(4)+d*y(2)*y(4)-yg80*y(4);
% dy(6) = a1*y(19)*y23-y1*y(6)-((vsd*y(6)*y(17))/(ksd+y(6)));
% dy(7) = a2*y22-y2*y(7);
% dy(8) = a3*y19*y21-y3*y(8)-((vsd*y(8)*y(17))/(ksd+y(8)));
% dy(9) = a80*y21-y80*y(9);
% dy(10)= k4d*(y(4)^2)-kn4d*y(10)-kr*y(11)*y(10)+knr*y(12)-d1*y(10);
% dy(11)= k80d*(y(5)^2)-kn80d*y(11)-kr*y(11)*y(10)+knr*y(12)-d1*y(11);
% dy(12)= kr*dy(11)*dy(10)-knr*y(12)-d1*y(12);
% dy(13)= k3*y(3)*y(15)-kn3*y(13)-k80*y(5)*y(13)+kn80*y(14)-d1*y(13);
% dy(14)= k80*y(5)*y(13)-kn80*y(14)-d1*y(14);
% dy(15)= y19*y20-((kcat*y(1)*y(15))/(kmgk+y(15)))-k3*y(3)*y(15)+kn3*y(13)-d1*y(15);
% dy(16)= (aglu+eglu*((y(18)/cglu)^b))/(1+((y(18)/cglu)^b))-yglu*y(16);
% dy(17)= tglu*y(16)-dglu*y(17);
% dy(18)= ktr2*y(17)*((GLUe-y(18))/(kmtr2+GLUe+y(18)+(atr2*GLUe*y(18))/kmtr2))-((uglu*y(18)*y(17))/(kglu+y(18)))-dd*y(18);




end



