clear all

file = 'M_20190408_04.txt';

%% ************************************************************************
%------------------------------GENERAL-------------------------------------
%**************************************************************************

% First get data's name
[~, name, ~] = fileparts(file);

% Now get data and its details
details = importdata(file, '', 7); % '.data' is data, '.textdata' is file's text
rawdata = importdata(file)
rawdata = rawdata.data;

% Get also some important parameters
nrepetitions = size(rawdata, 2) / 2; % Number of measurements
nsize = size(rawdata, 1);    % Length of each measurement

% Now create an equispaced time in ps
t = rawdata(:,1); % Out of simplicity
T = t(end) - t(1); % Total time
samplerate = nsize / T; % Sampling rate
dt = T / nsize; % Delta t
t = linspace(t(1), t(end), nsize)';

% And a matrix that holds the voltage in uV
for i = 1 : nrepetitions % 2 : nrepetitions in case of faulty measurement
    V(:,i) = 1e6 * rawdata(:, i*2);
end
meanV = mean(V, 2);

%% ************************************************************************
%-------------------------------PMUSIC-------------------------------------
%**************************************************************************

subplot(3,2,1);
    
nnn = floor(nsize / 4) % Number of measurements
Mp = [nnn, 200]; % Mean PMUSIC parameter
% Dimension, Limit that marks how many harmonics to discard
% Can't be greater than dimension
Mdt = 1200; % es el tamaño de la ventana y esta en ps
pasot = 20; % es el paso entre ventanas en ps

vmm = meanV;
vmm = detrend(vmm);
MSn1 = [];
Mfn1 = [];
MSn = [];
Mfn = [];
MC = [];
Mind = 0;

for i = Mdt+1 : pasot : 1350
    Mind = Mind+1;
    MSn1 = [];
    Mfn1 = [];
    seg = vmm(  ((i-Mdt)<t) & (t<i)  );
    [MSn1, Mfn1] = pmusic(seg, Mp, 6000, samplerate, [], 0);
    
    Sel = ( (Mfn1>=0.00) & (Mfn1<=0.06) );
    MSn(:,Mind) = MSn1(Sel);
    Mfn(:,Mind) = Mfn1(Sel);
end

MC = MSn;

imagesc([1:T], Mfn(:,1), MC);

%% ************************************************************************
%-------------------------LINEAR PREDICTION--------------------------------
%**************************************************************************

% Select data
x = meanV; % If I want to use the mean value.
%x = V(:,2); % If I want to use a particular measurement.

% Select t0
subplot(3,2,3);
plot(t, x);
xlim([-20, 1350]);
[t0, nn] = ginput;
t0

% Crop data
V = V(t>t0, :);
meanV = meanV(t>t0, :);
x = x(t>t0);
t = t(t>t0);

%%

% General parameters
c1=2e-5; c2=1e-5; c3=0.5e-4; cn=8e-15;
b1=3; b2=1; b3=5;
w1=5*2*pi; w2=4.5*2*pi;
p1=pi; p2=0.5*pi;

% Generate noise
for i=1:50;
    noise(i,:) = cos(1/(2*dt) * 2*pi * rand. * (t+dt)');
end
noise=cn*sum(noise);

% Coherent artifact (prop to croscorrelation)
coherent = 1e-15*exp(-t.^2 / 2 / 0.07^2);

N=length(t); 
%M=30;
M=round(0.75*N);

% Desired aspect of solution...
% x = c1.*exp(-b1.*t).*cos(w1*t+p1) + c2.*exp(-b2.*t).*cos(w2.*t+p2) + c3.*exp(-b3.*t.^2) + noise' + coherent;
% plot(t,x)

% Set up matrix from data (N-M)xM. We take M=0.75*N backward prediction.
X = zeros(N-M,M);
for i=1:M
    for k=1:N-M
        X(k,i)=x(i+k);
    end
end

% Computation of the (N-M)x(N-M) noonegativmatrix XX' and diagonalization
XX = X*X';
[U, D] = eig(XX);
d = eig(XX);

%%

no = 1:length(d);
%figure();
subplot(3,2,5);
semilogy(no,d,'o')
noo=input('no of significant Singular Values: ');
%close
%noo=8;

F=zeros(noo);

for j=length(d)-noo+1:length(d)
%for j=1:size(D,1)
  F(j,j)=1/sqrt(d(j));% 1 overlambda
    end


V=X'*U*F;


%computation of LP coeficients

xvector=x(1:N-M);
A=V*F'*U'*xvector;

%polynomial roots
%c=zeros(length(A));

  c(1)=1;

for i=1:length(A)
  c(i+1)=-A(i);
         end
         
         r=roots(c);
         

 %only l roots are significant being l rank of D matrix
         
         l=rank(D);

   s=sort(r);

    
 
   BB=length(s);
         
  %roots sorted in descending order

         for i=1:BB
 
          ss(i)=s(BB-i+1);
 
          end


%S=1;

   %for i=1:length(r)
   % if abs(ss(i))>=1
       %S=S+1;
     %end
     %end 
 
            ss=ss';
            
         for j=1:l

            b(j)=log(abs(ss(j)))/dt;
            w(j)=angle(ss(j))/dt;

                 end
              
         [P,I]=sort(w);
         Z=b(I);

       Nzeros=0;
         for j=1:l
           if w(j)==0
           Nzeros=Nzeros+1;
          end
         end
         Nzeros;

        for i=1:round((l-Nzeros)/2+Nzeros)
         WW(i)=abs(P(i));
         BBB(i)=(Z(i));  
         end
     %counting for positive damping constants  
     Npos=0;
       for j=1:length(BBB);
          if BBB(j)>=0
           Npos=Npos+1;
          end
        end

       [B1,J]=sort(BBB);
         W1=WW(J);


         
  % sorted in descending order

         for i=1:length(B1)
 
          B2(i)=B1(length(B1)-i+1);
 
          W2(i)=W1(length(B1)-i+1);

          end

        for j=1:Npos
           W(j)=W2(j);
          B(j)=B2(j);
          end    
         W;

         B;
         
         % LS for amplitudes and phases
         
         %setup matrices
   
        for i=1:N

            for j=1:length(W);
               Xbar(i,2*j-1)=exp(-B(j)*(i-1)*dt)*cos(W(j)*(i-1)*dt);
               Xbar(i,2*j)=-exp(-B(j)*(i-1)*dt)*sin(W(j)*(i-1)*dt);
             
            end
         end
        
     
           XXbar=Xbar*Xbar';
[Ubar,Dbar] = eig(XXbar);
dbar=eig(XXbar);

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



fi;
C;



for i=1:length(W)
   %if W(i)<40
      yy(:,i)=C(i).*exp(-B(i).*t).*cos(W(i).*t+fi(i));
   %else
   %      yy(:,i)=0*t;
   %end
end


yy=yy';
Y=sum(yy);
Chi2=sum((Y-x').^2)/N
%figure
%figure(2);
subplot(3,2,2);
plot(t,Y,t,x);title('Fit and Data');xlabel('time(ps)')


results=[W 
  B
  fi
  C];

results=results';

% up to here the programm, ahead calculations to evaluate if the fit is good enough


%calculation of the residue for evaluating the fit

Y=Y';
yy=yy';
residue=x-Y;
%figure;plot(t,residue);title('Residue');xlabel('time(ps)');ylabel('Data-Fit');

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

%spectrum of the fit in a "Ramanlike" fashion, meaning without considering the phase or phase zero

b=B/2/pi*100/3*30;%in GHz
f=W*100/3/2/pi*30;%in GHz
Wmax=max(W);
freq=0:Wmax/2/pi/1000*100/3*30:1.5*Wmax/2/pi*100/3*30;%in GHz
freq=freq';

res=zeros(length(freq),length(W));

for i=1:length(W)
   if W(i)==0
      res(:,i)=0;
   else 
      %res(:,i)=C(i)*imag(1i*exp(-1i*fi(i))*f(i)./(f(i)^2-freq.^2-2j.*freq*b(i)));
      res(:,i)=C(i)*imag(1*f(i)./(f(i)^2-freq.^2-2j.*freq*b(i)));
  end
  end
  
  spectrum=sum(res,2);
%figure;
subplot(3,2,6);
plot(freq,spectrum,freq,res);title('Raman-like Data Spectrum');xlabel('freq (GHz)');ylabel('spectrum');

%for SiO2
t02=0.26;
Cp=C.*exp(B*t02);
fip=W*t02-fi;

for i=1:length(W)
      yp(:,i)=Cp(i).*exp(-B(i).*t).*cos(W(i).*t+fip(i));
   end



%Los resultados en GHz
tau=1./B;%in picoseconds
fG=W*100/3/2/pi*30;%in GHz
Q=W/2./B;%factor de calidad

'Frecuencia en GHz, tau en ps, fase en grados, Amplitud, Q'
format shortE
resultsGHz=[fG 
   tau
   fi*180/pi
   C
   Q];
resultsGHz=resultsGHz'

%%

figure
plot(t,x)
hold on
plot(t,Y, 'LineWidth', 2.5)
xlabel('Tiempo (ps)')
ylabel('Voltaje (uV)')
savefig(data.name)

figure
% plot(freq,spectrum)
% hold on
plot(freq,res)
xlabel('Frecuencia (GHz)')
ylabel('Espectro reconstruído (u.a.)')
[~, len] = size(res);
for i=1:len
	legends(i) = sprintf('%f GHz', fG(i))
end
