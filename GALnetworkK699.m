
function dy = GALnetworkK699(t,y)
t
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
i=0;

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
GLUe = 0;%2.957*10^8;
p = 1.29*10^7;
q = 0.8;
vsd = 9.30*10^(-6);
ksd = 30;

y19= (p^q)/((p^q)+(y(17)^q));
y20= ktr*y(2)*((GALe-y(15))/(kmtr+GALe+y(15)+(atr*GALe*y(15)/kmtr)));

%% ----- For N=1 -------
N=1;
x=1;
a=1;
h=1;
i=0;
y(10)=1;
y(11)=1;
%Term of the nominator
for i= 0:N
  
     if (N==i)
           y2n(a) = 0;
           a =a+1;
     else
        for h=1:(N-i)

             y2n(a)= (nchoosek(N-i,h)*(cp^(h+i-1))*(cq^(i-1))*(kp*y(10))^h);
             a =a+1;
        end
     end 
 y1n(x)=(nchoosek(N,i)*(kq*kp*y(10)*y(11))^i);
 x=x+1;
end
yn= sum(y1n(1:2))*sum(y2n(1:2));
clear x; clear h; clear a;clear i;

%Term of the denominator
x=1;
a=1;

for i= 0:N
    
    for h=0:(N-i)
        y2d(a)= (nchoosek(1-i,h)*(cp^(h+i-1))*(cq^(i-1))*(kp*y(10))^h);
        a =a+1;
    end
 y1d(x)=(nchoosek(N,i)*(kq*kp*y(10)*y(11))^i);
 x=x+1;
end 
yd= sum(y1d(1:2))*sum(y2d(1:3));
clear x; clear h; clear a;
%Division
y21= yn/yd;

clear y2n;clear y1n;clear y2d; clear y1d;clear yd;clear yn;clear N;

%% ---- For N=2 -----
N=2;
x=1;
a=1;

%Term of the nominator
for i= 0:N
    
    if (N==i)
           y2n(a) = 0;
           a =a+1;
    else
       for h=1:(N-i)
          y2n(a)= (nchoosek(N-i,h)*(cp^(h+i-1))*(cq^(i-1))*(kp*y(10))^h);
          a =a+1;
       end
    end
 y1n(x)=(nchoosek(N,i)*(kq*kp*y(10)*y(11))^i);
 x=x+1;
end
yn= sum(y1n(1:3))*sum(y2n(1:3));
clear x; clear h; clear a;
%Term of the denominator
x=1;
a=1;
h=0;
for i= 0:N
    
   
    for h=0:(N-i)
        y2d(a)= (nchoosek(N-i,h)*(cp^(h+i-1))*(cq^(i-1))*(kp*y(10))^h);
        a =a+1;
    end
 y1d(x)=(nchoosek(N,i)*(kq*kp*y(10)*y(11))^i);
 x=x+1;
end
yd= sum(y1d(1:3))*sum(y2d(1:6));
clear x; clear h; clear a;
%Division
y22= yn/yd;

clear y2n;clear y1n;clear y2d; clear y1d;clear yd;clear yn;clear N;

%% ---- for N=4 -----

N=4;
x=1;
a=1;
h=1;
%Term of the nominator
for i= 0:N
   if (N==i)
           y2n(a) = 0;
           a =a+1;
   else
       for h=1:(N-i)
           y2n(a)= (nchoosek(N-i,h)*(cp^(h+i-1))*(cq^(i-1))*(kp*y(10))^h);
           a =a+1;
       end
   end 
 y1n(x)=(nchoosek(N,i)*(kq*kp*y(10)*y(11))^i);
 x=x+1;
end
yn= sum(y1n(1:5))*sum(y2n(1:5));
clear x; clear h; clear a;
%Term of the denominator
x=1;
a=1;
h=0;
for i= 0:N
   
    for h=0:(N-i)
        y2d(a)= (nchoosek(N-i,h)*(cp^(h+i-1))*(cq^(i-1))*(kp*y(10))^h);
        a =a+1;
    end
 y1d(x)=(nchoosek(N,i)*(kq*kp*y(10)*y(11))^i);
 x=x+1;
end
yd= sum(y1d(1:5))*sum(y2d(1:5));
clear x; clear h; clear a;
%Division
y23= yn/yd;

dy = zeros(23,1);

dy(1) = b1*y(6)-d1*y(1);
dy(2) = b2*y(7)-d1*y(2);
dy(3) = b3*y(8)-d1*y(2)-k3*y(15)*y(3)+kn3*y(13);
% the issue here is y19, is this y(19)?
dy(4) = b4*y19-d1*y(4)-2*k4d*((y(4))^2)+2*kn4d*y(10);
dy(5) = b80*y(9)-d1*y(5)-2*k80d*(y(5)^2)+2*kn80d*y(11)-k80*y(13)*y(5)+kn80*y(14);
dy(6) = a1*y(19)*y23-y1*y(6)-((vsd*y(6)*y(17))/(ksd+y(6)));
dy(7) = a2*y22-y2*y(7);
dy(8) = a3*y19*y21-y3*y(8)-((vsd*y(8)*y(17))/(ksd+y(8)));
dy(9) = a80*y21-y80*y(9);
dy(10)= k4d*(y(4)^2)-kn4d*y(10)-kr*y(11)*y(10)+knr*y(12)-d1*y(10);
dy(11)= k80d*(y(5)^2)-kn80d*y(11)-kr*y(11)*y(10)+knr*y(12)-d1*y(11);
dy(12)= kr*dy(11)*dy(10)-knr*y(12)-d1*y(12);
dy(13)= k3*y(3)*y(15)-kn3*y(13)-k80*y(5)*y(13)+kn80*y(14)-d1*y(13);
dy(14)= k80*y(5)*y(13)-kn80*y(14)-d1*y(14);
dy(15)= y19*y20-((kcat*y(1)*y(15))/(kmgk+y(15)))-k3*y(3)*y(15)+kn3*y(13)-d1*y(15);
dy(16)= (aglu+eglu*((y(18)/cglu)^b))/(1+((y(18)/cglu)^b))-yglu*y(16);
dy(17)= tglu*y(16)-dglu*y(17);
dy(18)= ktr2*y(17)*((GLUe-y(18))/(kmtr+GLUe+y(18)+(atr*GLUe*y(18))/kmtr))-((uglu*y(18)*y(17))/(kglu+y(18)))-dd*y(18);
%% FM 08/03/2016
% The ones below are not state variables but you are treating them as such.
% Let's discuss this tomorrow

% y19= (p^q)/((p^q)+(y(17)^q));
% y20= ktr*y(2)*((GALe-y(15))/(kmtr+GALe+y(15)+(atr*GALe*y(15)/kmtr)));
% y21= (nchoosek(N,i)*(kq*kp*y(10)*y(11))^i)*(nchoosek(1-i,h)*(cp^(h+i-1))*(cq^(i-1))*(kp*y(10))^h)/(nchoosek(1,i)*(kq*kp*y(10)*y(11))^i)*(nchoosek(1-i,h)*(cp^(h+i-1))*(cq^(i-1))*(kp*y(10))^h);

% Expressing the sum series with a for loop

%% This clear is not needed (none of these variables is exported from this function.
% clear y2n;clear y1n;clear y2d; clear y1d;clear yd;clear yn;



end



