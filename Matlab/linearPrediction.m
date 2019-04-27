clear all

rawdata = importdata('Datos.txt');
t = rawdata(1,:)';
x = rawdata(2,:)';
clear rawdata

%%

N = length(x);
M = round(N*0.75);

X = zeros(N-M,M);

for i=1:M
    for k=1:N-M
        X(k,i)=x(i+k);
    end
end
clear i k

%%

XX = X * X';
[U, D] = eig(XX);
d = eig(XX);

% for i = 1:length(d)
%     hey = U(1,i);
%     hello = U(:,i);
%     U(:,i) = hello/hey;
% end
% clear hello hey

%%

nsignificant = 4;

F = zeros(nsignificant);

for j=length(d)-nsignificant+1:length(d)
%for j=1:size(D,1)
  F(j,j)=1/sqrt(d(j));% 1 overlambda
end
    
%%

V=X'*U*F;
xvector=x(1:N-M);
A=V*F'*U'*xvector;

%%

c(1)=1;
for i=1:length(A)
    c(i+1)=-A(i);
end
r=roots(c);
s=sort(r); % sorted from smallest to largest absolute value
l=rank(D);

%%

for i=1:length(s)
    ss(i)=s(length(s)-i+1);    
end
ss = ss'; % sorted from largest to smallest absolute value

%%

dtime = 2;

for j=1:l
 
    b(j)=log(abs(ss(j)))/dtime;
    w(j)=angle(ss(j))/dtime;

end

%%

[P,I]=sort(w);
Z=b(I);

%%

Nzeros=0;
for j=1:l
    if w(j)==0
        Nzeros=Nzeros+1;
    end
end

for i=1:round((l-Nzeros)/2+Nzeros)
    WW(i)=abs(P(i));
    BBB(i)=(Z(i));  
end

%%

Npos=0;
for j=1:length(BBB);
    if BBB(j)>=0
        Npos=Npos+1;
    end
end

[B1,J]=sort(BBB);
W1=WW(J);

for i=1:length(B1)
    B2(i)=B1(length(B1)-i+1);
    W2(i)=W1(length(B1)-i+1);
end

for j=1:Npos
    W(j)=W2(j);
    B(j)=B2(j);
end   

%%

% LS for amplitudes and phases

% setup matrices

for i=1:N
    for j=1:length(W);
        Xbar(i,2*j-1)=exp(-B(j)*(i-1)*dtime)*cos(W(j)*(i-1)*dtime);
        Xbar(i,2*j)=-exp(-B(j)*(i-1)*dtime)*sin(W(j)*(i-1)*dtime);
    end
end


%%

XXbar = Xbar*Xbar';
[Ubar,Dbar] = eig(XXbar);
dbar = eig(XXbar);

%%

n1=1:length(dbar);
%figure
%plot(n1,dbar,'o')

%n11=input('number of points lambda');

n11=rank(XXbar);
Fbar=zeros(n11);
for j=length(dbar)-n11+1:length(dbar)

  Fbar(j,j)=1/sqrt(dbar(j));% 1 overlambda
    end


Vbar=Xbar'*Ubar*Fbar;


%least-squares

AA=Vbar*Fbar'*Ubar'*x;

%%

for i=1:length(W)

     if AA(2*i-1)==0 & AA(2*i)==0
           C(i)=0;fi(i)=0;
     elseif AA(2*i-1)==0
            fi(i)=sign(AA(2*i))*pi/2;
             C(i)=abs(AA(2*i));
     elseif AA(2*i)==0
             fi(i)=(1-sign(AA(2*i-1)))*pi/2;
             C(i)=abs(AA(2*i-1));
        else
   fi(i)=atan2(AA(2*i),AA(2*i-1));
 
    C(i)=sqrt(AA(2*i)^2+AA(2*i-1)^2);
    end
    end

%%

for i=1:length(W)
   %if W(i)<40
      yy(:,i)=C(i).*exp(-B(i).*(t-t(1))).*cos(W(i).*(t-t(1))+fi(i));
   %else
   %      yy(:,i)=0*t;
   %end
end

yy = yy';
Y = sum(yy);
Chi2=sum((Y-x').^2)/N

%%

% up to here the programm, ahead calculations to evaluate if the fit is good enough

%calculation of the residue for evaluating the fit

residue=x'-Y;

%FFt of the residue

%Y1 = fft(residue,4096);
Y1=fft(residue);
NN=length(Y1);
Pyy= Y1.*conj(Y1)/(NN-1);
frequ = 1/(t(2)-t(1))/NN*(0:((NN/2)-1))*100/3*30;%in GHz
FFTResidue=Pyy(1:NN/2);
frequ=frequ';
%figure;
subplot(3,2,4);
plot(frequ,Pyy(1:NN/2));title('FFT of the Residue');xlabel('freq (GHz)');ylabel('FFT');
