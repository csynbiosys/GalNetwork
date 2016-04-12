% [t,y]=ode45(@GALnetworkK699,[0 1000],ones(7,1))
function dy = GALnetworkK699(t,y)
t


% new values from the simplified model
s = 100 ;%2.366*10^8; %Galactose concentration outside the cell
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




dy = zeros(7,1);
%Equations of the glucose module
dy(5)= ktr2*y(7)*((GLUe-y(5))/(kmtr2+GLUe+y(5)+(atr2*GLUe*y(5))/kmtr2))-((uglu*y(5)*y(7))/(kglu+y(5)))-dd*y(5); %Internal glucose concentration
dy(6)= (aglu+eglu*((y(5)/cglu)^b))/(1+((y(5)/cglu)^b))-yglu*y(6); %Glucose network mRNA
dy(7)= tglu*y(6)-dglu*y(7);%Glucose network protein
y19= (p^q)/((p^q)+(y(7)^q)); %This is the glucose repression function


% Equations for the protein production
dy(1) = a0g1+as*e+ag1*((y(3)^n1)/((y(3)^n1)+Kg1^n1))+w*y(1)*y(4)-yg1*y(1); %Protein concentration of Gal1
dy(2) = as + ag3*((y(3)^n3)/((y(3)^n3)+Kg3^n3))+d*y(2)*y(4)-yg3*y(2); %Protein concentration of Gal3
dy(3) = ag4*y19+b*y(3)*y(4)-yg4*y(3); %Protein concentration of Gal4
dy(4) = a0g80+ag80*((y(3)^n80)/((y(3)^n80)+Kg80^n80))+w*y(1)*y(4)+b*y(3)*y(4)+d*y(2)*y(4)-yg80*y(4);%Protein concentration of Gal80






end



