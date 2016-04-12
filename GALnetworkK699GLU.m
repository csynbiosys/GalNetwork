% [t,y]=ode45(@GALnetworkK699GLU,[0 1000],ones(3,1))
function dy = GALnetworkK699GLU(t,y)
t
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
GLUe = 2.957*10^8;
p = 1.29*10^7;
q = 0.8;

dy = zeros(3,1);
%Equations of the glucose module
dy(1)= ktr2*y(3)*((GLUe-y(1))/(kmtr2+GLUe+y(1)+(atr2*GLUe*y(1))/kmtr2))-((uglu*y(1)*y(3))/(kglu+y(1)))-dd*y(1);
dy(2)= (aglu+eglu*((y(1)/cglu)^b))/(1+((y(1)/cglu)^b))-yglu*y(2);
dy(3)= tglu*y(2)-dglu*y(3);
y19= (p^q)/((p^q)+(y(3)^q));

end




