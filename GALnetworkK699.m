
function dy = GALnetworkK699(t,y)
% values for galactose network module for the K699-strain
b1 = 9.92;
b2 = 6.94;
b3 = 18.0;
b4 = 0.86;
b80= 4.00;
a1 = 1.09;
a2 = 1.20;
a3 = 36.0;
a80= 3.00;
d1 = 0.0033;
y1 = 0.036;
y2 = 0.026;
y3 = 0.036;
y80 = 0.036;
k3 = 5.0*10^(-8);
kn3 = 890;
k4d = 0.10;
kn4d = 1.0;
k80d = 0.10;
kn80d = 170;
k80 = 0.10;
kn80 = 0.03;
kr = 0.10;
knr = 1.80;
kcat = 3350;
kmgk = 1.29*10^7;
ktr = 4350;
kmtr = 2.15*10^8;
atr = 10.0;
kp = 0.091;
kq = 0.0556;
cp = 1;
cq = 30;
GALe = 2.366*10^8;

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
kmtr = 6.022*10^8;
atr = 1.0;
uglu = 5350;
kglu = 1.29*10^7;
GLUe = 2.957*10^8;
p = 1.29*10^7;
q = 0.8;
vsd = 9.30*10^(-6);
ksd = 30;


dy = zeros(18,1);
dy(1) = b1*y(6)-d1*y(1);
dy(2) = b2*y(7)-d1*y(2);
dy(3) = b3*y(8)-d1*y(2)-k3*y(15)*y(3)+kn3*y(13);
dy(4) = b4*y(19)-d1*y(4)-2*k4d*((y(4))^2)+2*kn4d*y(10);
dy(5) = b80*y(9)-d1*y(5)-2*k80d*(y(5)^2)+2*kn80d*y(11)-k80*y(13)*y(5)+kn80*y(14);
dy(6) = a1*y(19)*y(23)-y1*y(6)-((vsd*y(6)*y(17))/(ksd+y(6)));
dy(7) = a2*y(22)-y2*y(7);
dy(8) = a3*y(19)*y(21)-y3*y(8)-((vsd*y(8)*y(17))/(ksd+y(8)));
dy(9) = a80*y(21)-y80*y(9);
dy(10)= k4d*(y(4)^2)-kn4d*y(10)-kr*y(11)*y(10)+knr*y(12)-d1*y(10);
dy(11)= k80d*(y(5)^2)-kn80d*y(11)-kr*y(11)*y(10)+knr*y(12)-d1*y(11);
dy(12)= kr*dy(11)*dy(10)-knr*y(12)-d1*y(12);
dy(13)= k3*y(3)*y(15)-kn3*y(13)-k80*y(5)*y(13)+kn80*y(14)-d1*y(13);
dy(14)= k80*y(5)*y(13)-kn80*y(14)-d1*y(14);
dy(15)= y(19)*y(20)-((kcat*y(1)*y(15))/(kmgk+y(15)))-k3*y(3)*y(15)+kn3*y(13)-d1*y(15);
dy(16)= (aglu+eglu*((y(18)/cglu)^b))/(1+((y(18)/cglu)^b))-yglu*y(16);
dy(17)= tglu*y(16)-dglu*y(17);
dy(18)= ktr2*y(17)*((GLUe-y(18))/(kmtr2+GLUe+y(18)+(atr2*GLUe*y(18))/kmtr2))-((uglu*y(18)*y(17))/(kglu+y(18)))-dd*y(18);

y(19)= (p^q)/((p^q)+(y(17)^q));
y(20)= ktr*y(2)*((GALe-y(15))/(kmtr+GALe+y(15)+(atr*GALe*y(15)/kmtr)));
% For N=1
y(21)= (symsum((nchoosek(1,i)*(kq*kp*y(10)*y(11))^i),i,0,1))*(symsum((nchoosek(1-i,h)*(cp^(h+i-1))*(cq^(i-1))*(kp*y(10))^h),h,1,1-i))/((symsum((nchoosek(1,i)*(kq*kp*y(10)*y(11))^i),i,0,1))*(symsum((nchoosek(1-i,h)*(cp^(h+i-1))*(cq^(i-1))*(kp*y(10))^h),h,0,1-i)));
% For N=2
y(22)= (symsum((nchoosek(2,i)*(kq*kp*y(10)*y(11))^i),i,0,2))*(symsum((nchoosek(2-i,h)*(cp^(h+i-1))*(cq^(i-1))*(kp*y(10))^h),h,1,2-i))/((symsum((nchoosek(2,i)*(kq*kp*y(10)*y(11))^i),i,0,2))*(symsum((nchoosek(2-i,h)*(cp^(h+i-1))*(cq^(i-1))*(kp*y(10))^h),h,0,2-i)));
% For N=4
y(23)= (symsum((nchoosek(4,i)*(kq*kp*y(10)*y(11))^i),i,0,4))*(symsum((nchoosek(4-i,h)*(cp^(h+i-1))*(cq^(i-1))*(kp*y(10))^h),h,1,4-i))/((symsum((nchoosek(4,i)*(kq*kp*y(10)*y(11))^i),i,0,4))*(symsum((nchoosek(4-i,h)*(cp^(h+i-1))*(cq^(i-1))*(kp*y(10))^h),h,0,4-i)));


end



