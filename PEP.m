clc;
r=1:0.1:10;
l=length(r);

e=zeros(1,l);

Mt=2;
Mr=1;
C0=[1 1;-1 1];
C1=[-1 -1;1 -1];

c=sqrt(rdivide(r,Mt));
N=1000000;
for n=1:N
    
h1=randn(1)+1i*randn(1);
h2=randn(1)+1i*randn(1);
H=[h1;h2];

n01=randn(1)+1i*randn(1);
n02=randn(1)+1i*randn(1);
N0=[n01;n02];

n1=randn(1)+1i*randn(1);
n2=randn(1)+1i*randn(1);
N1=[n1;n2];

C0H= (C0*H);
C1H= (C1*H);
Y0=(c.*C0H) +N0;
Y1=(c.*C1H) +N1;
s=[(Y0-(c.*C0H)); (Y0-(c.*C1H))];

n1=zeros(1,l);
n2=zeros(1,l);
i=0;

for i=0:l-1
n1(i+1)=norm(s(1:2,i+1),'fro');
end

for i=0:l-1
n2(i+1)=norm(s(3:4,i+1),'fro');
end
p1=n1.^2;
p2=n2.^2;

p3=p1-p2;


for i=0:l-1
    if (p3(i+1)>0)
        e(i+1)=e(i+1)+1;
    end
end
end
e;
pe=rdivide(e,N);
lpe= 10*log10(pe);

x=sqrt(2*r*(((abs(h1))^2) + ((abs(h2))^2)));
ipep=qfunc(x);
tpep=0.5* ((rdivide(1,(1+r))).^2);
tlpep= 10*log10(tpep);
lr=10*log10(r);
%subplot(1,2,1)
semilogy(lr,tpep,lr,pe);

%subplot(1,2,2)
%semilogy(r,lpe);