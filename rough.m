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
h6(:,l)=20 * (log10(h8(:,l)));

end


%{
subplot(3,3,1);
semilogy(k,h6(:,1));
grid on;
title('ques 2');
xlabel('Channel 1 samples [K]');
ylabel('Amplitude in dB');

subplot(3,3,2);
semilogy(k,h6(:,2));
grid on;
title('ques 2');
xlabel('Channel 2 samples [K]');
ylabel('Amplitude in dB');

subplot(3,3,3);
semilogy(k,h6(:,3));
grid on;
title('ques 2');
xlabel('Channel 3 samples [K]');
ylabel('Amplitude in dB');

subplot(3,3,4);
semilogy(k,h6(:,4));
grid on;
title('ques 2');
xlabel('Channel 4 samples [K]');
ylabel('Amplitude in dB');


subplot (3,3,5);
N1=histcounts(h8(:,1),'Normalization','pdf');
plot(N1);
grid on;
xlabel('Envelope of the Experimental pdf for 2nd ques');

subplot (3,3,6);
N1=histcounts(h8(:,2),'Normalization','pdf');
plot(N1);
grid on;
xlabel('Envelope of the Experimental pdf for 2nd ques');

subplot (3,3,7);
N1=histcounts(h8(:,3),'Normalization','pdf');
plot(N1);
grid on;
xlabel('Envelope of the Experimental pdf for 2nd ques');

subplot (3,3,8);
N1=histcounts(h8(:,4),'Normalization','pdf');
plot(N1);
grid on;
xlabel('Envelope of the Experimental pdf for 2nd ques');


subplot(3,3,9);
%}
k1=1:1:500;
plot(k1,h6(1:500,1),':r',k1,h6(1:500,2),':b*',k1,h6(1:500,3),'--g',k1,h6(1:500,4),'black');
grid on;
title('ques 2');
xlabel('Channel samples [K]');
ylabel('Amplitude |h(t)| in dB');
