clc;
clear all; % a must 
%to generate 4 independent rayleigh fading channels
fdm= 59.60518;
Ts=1/(3*(10^(4)));
fdmn= fdm * Ts;

M=8;
N=34;
K=100000;
k=1:1:K;
n=1:1:M;
L=4;

h= zeros(K,L);
h1= zeros(K,L);
h2= zeros(K,L);
h3= zeros(K,L);
h4= zeros(K,L);
h5= zeros(K,L);
h6= zeros(K,L);
h8= zeros(K,L);
E0= zeros(K,L);

had8= hadamard(8);

Bn= (pi * n) / (M+1);
thn = ( 2 * pi * n ) / N;

for l=1:1:L
Gnl= (2 * pi * (l+1)*n )/( M+1);
Gnl1= padarray(Gnl,[K-1 0],'replicate','post');

c= (2 * pi * (k.') * fdmn * cos(thn))+ Gnl1;

y=cos(c);
e1=(exp(1i*Bn));
e = e1.*had8(l,:);

p= zeros(K,M);
for o=1:K
p(o,:) = y(o,:) .*e;
end

h(:,l) = sum(p,2);


%h1(:,l) = 10 * log10(h(:,l));
h2(:,l) = abs (h(:,l));
h3(:,l) = (h2(:,l)).^2;
h4(:,l) = sum(h3(:,l));
E0(:,l) = h4(:,l)/K;
h5(:,l) = h(:,l) ./ (sqrt(E0(:,l)));
h8(:,l) = abs(h5(:,l));
h6(:,l)= 10 * (log10(h8(:,l)));

end


%2(ii)
Mt=2;
Imt= eye(Mt);
S0= (sqrt(Mt)) * Imt;
Cl1= sqrt(2) * Imt * (exp((1i*0*pi)/2));
Cl2= sqrt(2) * Imt * (exp((1i*1*pi)/2));
Cl3= sqrt(2) * Imt * (exp((1i*2*pi)/2));
Cl4= sqrt(2) * Imt * (exp((1i*3*pi)/2));
SER = zeros(10,K-1);
SER_1 = zeros(10,K-1);
for repeat= 1:1:1
    
Nt=(1/sqrt(2)).*[(randn(1)+1i*randn(1)) ;(randn(1)+1i*randn(1)) ];
Nt1=(1/sqrt(2)).*[(randn(1)+1i*randn(1)) (randn(1)+1i*randn(1));(randn(1)+1i*randn(1)) (randn(1)+1i*randn(1))];

j=1;
for r=1:1:10
St_1 = S0;
Ht=[h5(1,1);h5(1,2)];
Ht1=[h5(1,1) h5(1,3);h5(1,2) h5(1,4)];
St= rdivide((Cl1 * St_1),sqrt(Mt));
Yt= rdivide((sqrt(r) * (St* Ht)),sqrt(Mt)) + Nt;
Yt1= rdivide((sqrt(r) * (St* Ht1)),sqrt(Mt)) + Nt1;
Yt1_1=Yt1;
Yt_1=Yt;
St_1=St;
Nt=(1/sqrt(2)).*[(randn(1)+1i*randn(1)) ;(randn(1)+1i*randn(1)) ];
Nt1=(1/sqrt(2)).*[(randn(1)+1i*randn(1)) (randn(1)+1i*randn(1));(randn(1)+1i*randn(1)) (randn(1)+1i*randn(1))];


i=1;
for k=2:1:K
Ht=[h5(k,1);h5(k,2)];
Ht1=[h5(k,1) h5(k,3);h5(k,2) h5(k,4)];
St= rdivide((Cl1 * St_1),sqrt(Mt));
Yt = rdivide((sqrt(r) * (St * Ht)),sqrt(Mt)) + Nt;
Yt1= rdivide((sqrt(r) * (St* Ht1)),sqrt(Mt)) + Nt1;
clt(1)=(norm((Yt - rdivide((Cl1*Yt_1),sqrt(Mt))),'fro'))^2;
clt(2)=(norm((Yt - rdivide((Cl2*Yt_1),sqrt(Mt))),'fro'))^2;
clt(3)=(norm((Yt - rdivide((Cl3*Yt_1),sqrt(Mt))),'fro'))^2;
clt(4)=(norm((Yt - rdivide((Cl4*Yt_1),sqrt(Mt))),'fro'))^2;

clt1(1)=(norm((Yt1 - rdivide((Cl1*Yt1_1),sqrt(Mt))),'fro'))^2;
clt1(2)=(norm((Yt1 - rdivide((Cl2*Yt1_1),sqrt(Mt))),'fro'))^2;
clt1(3)=(norm((Yt1 - rdivide((Cl3*Yt1_1),sqrt(Mt))),'fro'))^2;
clt1(4)=(norm((Yt1 - rdivide((Cl4*Yt1_1),sqrt(Mt))),'fro'))^2;

[M,I]=min(clt(:));
[M1,I1]=min(clt1(:));
if(I~=1)
    SER(j,i) = SER(j,i) +1;
end

if(I1~=1)
    SER_1(j,i) = SER_1(j,i) +1;
end

Nt=(1/sqrt(2)).*[(randn(1)+1i*randn(1)) ;(randn(1)+1i*randn(1)) ];
Nt1=(1/sqrt(2)).*[(randn(1)+1i*randn(1)) (randn(1)+1i*randn(1));(randn(1)+1i*randn(1)) (randn(1)+1i*randn(1))];

Yt_1=Yt;
Yt1_1=Yt1;
St_1=St;
i=i+1;
end
j=j+1;
end
end

SER1= rdivide(SER,1);
MSER= mean(SER1,2)

SER1_1= rdivide(SER_1,1);
MSER1= mean(SER1_1,2)

r=1:1:10;
SNR= 10*log10(r);

semilogy(SNR,MSER,':',SNR,MSER1,':r*');

