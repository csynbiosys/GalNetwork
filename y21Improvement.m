N=1;
x=1;
a=1;
h=1;
i=0;
y(10)=1;
y(11)=1;
cp =1;
cq=30;
kp = 0.091;
kq= 0.0556;
%Term of the nominator
% for i= 0:N
%   
%      if (N==i)
%            y2n(a) = 0;
%            a =a+1;
%      else
%         for h=1:(N-i)
% 
%              y2n(a)= (nchoosek(N-i,h)*(cp^(h+i-1))*(cq^(i-1))*(kp*y(10))^h);
%              a =a+1;
%         end
%      end 
%  y1n(x)=(nchoosek(N,i)*(kq*kp*y(10)*y(11))^i);
%  x=x+1;
% end

i = 0:1;
ylength = length(i);
y1n(1:ylength)=(nchoosek(1,i)*(kq*kp*y(10)*y(11))^i);

% yn= sum(y1n(1:2))*sum(y2n(1:2));
 clear x; clear h; clear a;clear i;

% %Term of the denominator
% x=1;
% a=1;
% 
% for i= 0:N
%     
%     for h=0:(N-i)
%         y2d(a)= (nchoosek(1-i,h)*(cp^(h+i-1))*(cq^(i-1))*(kp*y(10))^h);
%         a =a+1;
%     end
%  y1d(x)=(nchoosek(N,i)*(kq*kp*y(10)*y(11))^i);
%  x=x+1;
% end 
% yd= sum(y1d(1:2))*sum(y2d(1:3));
% clear x; clear h; clear a;
% %Division
% y21= yn/yd;
