clc;
clear all;
L=32;
l=0:1:(L-1);
l2=1:1:L-1;
Mt=2;
Thl=(l.*(2*pi))/L;
Thl2=(l2.*pi)/L;
i=1;
u11=0:0.01:L-1;
u12=0:0.01:L-1;
u=zeros(2,4);

for u1=0:0.01:L-1

    
for u2=0:0.01:L-1
a=((abs((sin(u1.*Thl2)).*(sin(u2.*Thl2)) )).^(1/Mt));
m=min(a);
ma(1,i)=m;
i=i+1;
end


end
L1=length(u11);
L2=length(u12);

[M,I]=max(ma)
y=I/L1;
x=mod(I,L1);
y1=y-(x/L1);
if x==0
u02=(L1-1)*u11(1,2)
else
    u02=(x-1)*u11(1,2)
end
u01=y1*u11(1,2)







%sin(rdivide((U1'.*Thl2),2));
%sin(rdivide((U2'.*Thl2),2));
