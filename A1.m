clc;
fdm= 80.467;
Ts=1/(3*(10^(4)));
fdmn= fdm * Ts

M=8;
N=34;
K=100000;
k=1:1:K;
n=1:1:M;
Bn= (pi * n) / (M+1);
thn = ( 2 * pi * n ) / N;
c= 2 * pi * (k.') * fdmn * cos(thn);

y=cos(c);

e = exp(1i*Bn);

p= zeros(K,M);
for o=1:K;
p(o,:) = y(o,:) .*e;
end;

h = sum(p,2);


h1 = 10 * log10(h);
h2 = abs (h1);
h3 = (h2).^2;
h4 = sum(h3(:));
E0 = h4/K;
h5 = h ./ (sqrt(E0));
h8 = abs(h5);
h6= 10 * (log10(h8));
%h7 = abs(h6); 


fdmn1=0.1;
c1= 2 * pi * (k.') * fdmn1 * cos(thn);
y1=cos(c1);
p1 = y1 .*e;
h2_0= sum(p1,2);
h2_1 = 10 * log10(h2_0);
h2_2 = abs (h2_1);
h2_3 = (h2_2).^2;
h2_4 = sum(h2_3(:));
E2_0 = h2_4/K;
h2_5 = h2_0 ./ (sqrt(E2_0));
h2_8 = abs(h2_5);
h2_6 = 10 * (log10(h2_8));
%h2_7 = abs(h2_6);


subplot(3,3,7);
semilogy(k,h2_6);
grid on;
title('ques 3');
xlabel('Channel samples [K]');
ylabel('Amplitude in dB');

subplot(3,3,6);
semilogy(k,h6);
grid on;
title('ques 2');
xlabel('Channel samples [K]');
ylabel('Amplitude in dB');

x7= 0:0.001:3;
y7 = raylpdf(x7,1/sqrt(2));
subplot(3,3,1);
plot(x7,y7);
grid on;
xlabel('Theoretical rayleigh pdf');

subplot(3,3,2);
histogram(h8,200,'Normalization','pdf');
grid on;
xlabel('Experimental pdf for 2nd ques');

subplot (3,3,4);
N1=histcounts(h8,'Normalization','pdf');
plot(N1);
grid on;
xlabel('Envelope of the Experimental pdf for 2nd ques');

subplot (3,3,5);
N2_1=histcounts(h2_8,'Normalization','pdf');
plot(N2_1);
grid on;
xlabel('Envelope of the Experimental pdf Fdmax.Ts =0.1');

subplot (3,3,3);
histogram(h2_8,200,'Normalization','pdf');
grid on;
xlabel(' Experimental pdf for Fdmax.Ts =0.1');
