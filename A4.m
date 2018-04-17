clc;

r=1:0.1:12;
l=length(r);
repeat=1:1:10000;
len1=length(repeat);
Mt=2;
x1=[1 -1 1i -1i];
x2=[1 -1 1i -1i];

l1=length(x1);
l2=length(x2);
e=zeros(l,1);
o=1;
for repeat=1:1:10000


n1=randn(1)+1i*randn(1);
n2=randn(1)+1i*randn(1);
N=[n1;n2];

h1=randn(1)+1i*randn(1);
h2=randn(1)+1i*randn(1);
H=[h1;h2];



k=1;
for r=1:0.1:12 
  j=1;  
 L1= 1;
 L2= 1;

G2=[x1(L1),x2(L2);-conj(x2(L2)),conj(x1(L1))];
Y=((sqrt(r/Mt)).*(G2*H)) + N;
i=1;
d1=zeros(1,16);
for L11=1:1:length(x1)
    for L12=1:1:length(x2)
        G2=[x1(L11),x2(L12);-conj(x2(L12)),conj(x1(L11))];
        d= Y-((sqrt(r/Mt)).*(G2*H));
        d1(1,i)=(norm(d,'fro'))^2;
        i=i+1;
        
    end
end
d1;
[M,I]=min(d1);

y=I/l1;
x=mod(I,l1);
y1=y-(x/l1);
if x==0
L01=y;
L02=l1;
else
 L02=x;
 L01=y1+1;
end

if((L1~=L01)||(L2~=L02))
    if((L1~=L01))
        if(((L1==1)&&(L01==4))||((L01==1)&&(L1==4)))
        e(k,j)=e(k,j)+2;
        else
            if(((L1==2)&&(L01==3))||((L01==2)&&(L1==3)))
        e(k,j)=e(k,j)+2;
            else 
              e(k,j)=e(k,j)+1;  
            end
            
        end
    end
    
       if((L2~=L02))
        if(((L2==1)&&(L02==4))||((L02==1)&&(L2==4)))
        e(k,j)=e(k,j)+2;
        else
            if(((L2==2)&&(L02==3))||((L02==2)&&(L2==3)))
        e(k,j)=e(k,j)+2;
            else 
              e(k,j)=e(k,j)+1;  
             end
        end
        
       end
       
end

k=k+1;
end
end
be=rdivide(e,4*len1)
be2 =( mean(be,2))
r=1:0.1:12;
r2b=1:0.1:12;
l2b=length(r2b);
repeat2b=1:1:10000;
len12b=length(repeat2b);
x12b=[1 -1 1i -1i];

l12b=length(x12b);

e2b=zeros(l2b,1);

for repeat2b=1:1:10000

n12b=randn(1)+1i*randn(1);

h12b=randn(1)+1i*randn(1);

j2b=1;
for r2b=1:0.1:12
 
L12b=1;
 
Y2b=((sqrt(r2b))*(x12b(L12b)*h12b)) + n12b;
i2b=1;
d12b=zeros(1,4);
for L11=1:1:length(x12b)
    
        d2b= Y2b-((sqrt(r2b))*(x12b(L11)*h12b));
        d12b(1,i2b)=(norm(d2b,'fro'))^2;
        i2b=i2b+1;
end
[M2b,I2b]=min(d12b);
L012b=I2b;

if((L12b~=L012b))
        if(((L12b==1)&&(L012b==4))||((L012b==1)&&(L12b==4)))
        e2b(j2b,1)=e2b(j2b,1)+2;
        else
            if(((L12b==2)&&(L012b==3))||((L012b==2)&&(L12b==3)))
        e2b(j2b,1)=e2b(j2b,1)+2;
            else 
              e2b(j2b,1)=e2b(j2b,1)+1;
            end
        end   
end
j2b=j2b+1;
end

end
be2b=rdivide(e2b,2*len12b);
be22b=(mean(be2b,2));
r2b=1:0.1:12;




semilogy(r,be2,'r',r2b,be22b,'g');

