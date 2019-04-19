# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:28:53 2019

@author: LEC
"""

import numpy as np
from numpy import pi
import matplotlib.pyplot as plt

#%%

def linearPrediction(x, t, dt):

    # Define general parameters
    c1 = 2e-5
    c2 = 1e-5
    c3 = 0.5e-4
    cn = 8e-15
     
    b1 = 3
    b2 = 1
    b3 = 5
    
    w1 = 5*2*pi
    w2 = 4.5*2*pi
    
    p1 = pi
    p2 = 0.5*pi
    
    # Generate noise and coherent artifect
    noise = []
    for i in range(50):
        noise.append(np.cos(2*pi*np.random.rand() * (t+dt) / (2*dt)))
    noise = np.array(noise)
    noise = cn * sum(noise)
    coherent_artifact = 1e-15 * np.exp(-t**2 / 2 / 0.07**2)
    # Proporcional to croscorrelation
    
    """
    Kind of answer we want:
    x = c1.*exp(-b1.*t).*cos(w1*t+p1) + c2.*exp(-b2.*t).*cos(w2.*t+p2) + 
        + c3.*exp(-b3.*t.^2) + noise + coherent_artifact;
    """

    # Set up matrix from data (N-M) x M. We take M=0.75*N backward prediction
    N = len(t)
    M = round(0.75 * N)
    X = np.zeros( (N-M, M) )
    for i in range(M):
        for k in range(N-M):
            X[k,i] = x[i+k+1]

    # Computation of the (N-M)x(N-M) noonegativmatrix XX' and diagonalization
    [ev, eV] = np.linalg.eig( np.matmul(X, X.T) )
    ordered_index = ev.argsort()
    ev = ev[ordered_index] # eigenvalues
    eV = eV[:, ordered_index] # eigenvectors on columns
    eD = np.diag(ev) # diagonal containing eigenvalues

    # XX is np.matmul(X, X.T)
    # U is ev eigenvalues
    # D is eV eigenvectors
    # noo is nsignificant

    # Choose number of singular values    
    plt.figure()
    plt.semilogy(ev, linestyle='none', marker='o', 
                 fillstyle='none', markersize=10)
    plt.xlabel("Número")
    plt.ylabel("Autovalores")
    nsignificant = int(input('¿Número de valores singulares?\t'))
    
    F = np.diag(1/np.sqrt(ev[:nsignificant]))
    U = np.matmul(X, eV*F)

#    #computation of LP coeficients
#    
#    xvector=x(1:N-M);
#    A=V*F'*U'*xvector;
#    
#    #polynomial roots
#    #c=zeros(length(A));
#    
#      c(1)=1;
#    
#    for i=1:length(A)
#      c(i+1)=-A(i);
#             end
#             
#             r=roots(c);
#             
#    
#     #only l roots are significant being l rank of D matrix
#             
#             l=rank(D);
#    
#       s=sort(r);
#    
#        
#     
#       BB=length(s);
#             
#      #roots sorted in descending order
#    
#             for i=1:BB
#     
#              ss(i)=s(BB-i+1);
#     
#              end
#    
#    
#    #S=1;
#    
#       #for i=1:length(r)
#       # if abs(ss(i))>=1
#           #S=S+1;
#         #end
#         #end 
#     
#                ss=ss';
#                
#             for j=1:l
#    
#                b(j)=log(abs(ss(j)))/dtime;
#                w(j)=angle(ss(j))/dtime;
#    
#                     end
#                  
#             [P,I]=sort(w);
#             Z=b(I);
#    
#           Nzeros=0;
#             for j=1:l
#               if w(j)==0
#               Nzeros=Nzeros+1;
#              end
#             end
#             Nzeros;
#    
#            for i=1:round((l-Nzeros)/2+Nzeros)
#             WW(i)=abs(P(i));
#             BBB(i)=(Z(i));  
#             end
#         #counting for positive damping constants  
#         Npos=0;
#           for j=1:length(BBB);
#              if BBB(j)>=0
#               Npos=Npos+1;
#              end
#            end
#    
#           [B1,J]=sort(BBB);
#             W1=WW(J);
#    
#    
#             
#      # sorted in descending order
#    
#             for i=1:length(B1)
#     
#              B2(i)=B1(length(B1)-i+1);
#     
#              W2(i)=W1(length(B1)-i+1);
#    
#              end
#    
#            for j=1:Npos
#               W(j)=W2(j);
#              B(j)=B2(j);
#              end    
#             W;
#    
#             B;
#             
#             # LS for amplitudes and phases
#             
#             #setup matrices
#       
#            for i=1:N
#    
#                for j=1:length(W);
#                   Xbar(i,2*j-1)=exp(-B(j)*(i-1)*dtime)*cos(W(j)*(i-1)*dtime);
#                   Xbar(i,2*j)=-exp(-B(j)*(i-1)*dtime)*sin(W(j)*(i-1)*dtime);
#                 
#                end
#             end
#            
#         
#               XXbar=Xbar*Xbar';
#    [Ubar,Dbar] = eig(XXbar);
#    dbar=eig(XXbar);
#    
#    n1=1:length(dbar);
#    #figure
#    #plot(n1,dbar,'o')
#    
#    #n11=input('number of points lambda');
#    
#    n11=rank(XXbar);
#    Fbar=zeros(n11);
#    for j=length(dbar)-n11+1:length(dbar)
#    
#      Fbar(j,j)=1/sqrt(dbar(j));# 1 overlambda
#        end
#    
#    
#    Vbar=Xbar'*Ubar*Fbar;
#    
#    
#    #least-squares
#    
#    AA=Vbar*Fbar'*Ubar'*x;
#    
#    
#    
#    for i=1:length(W)
#    
#         if AA(2*i-1)==0 & AA(2*i)==0
#               C(i)=0;fi(i)=0;
#         elseif AA(2*i-1)==0
#                fi(i)=sign(AA(2*i))*pi/2;
#                 C(i)=abs(AA(2*i));
#         elseif AA(2*i)==0
#                 fi(i)=(1-sign(AA(2*i-1)))*pi/2;
#                 C(i)=abs(AA(2*i-1));
#            else
#       fi(i)=atan2(AA(2*i),AA(2*i-1));
#     
#        C(i)=sqrt(AA(2*i)^2+AA(2*i-1)^2);
#        end
#        end
#    
#    
#    
#    fi;
#    C;
#    
#    
#    
#    for i=1:length(W)
#       #if W(i)<40
#          yy(:,i)=C(i).*exp(-B(i).*t).*cos(W(i).*t+fi(i));
#       #else
#       #      yy(:,i)=0*t;
#       #end
#    end
#    
#    
#    yy=yy';
#    Y=sum(yy);
#    Chi2=sum((Y-x').^2)/N
#    #figure
#    #figure(2);
#    subplot(3,2,2);
#    plot(t,Y,t,x);title('Fit and Data');xlabel('time(ps)')
#    
#    
#    results=[W 
#      B
#      fi
#      C];
#    
#    results=results';
#    
#    # up to here the programm, ahead calculations to evaluate if the fit is good enough
#    
#    
#    #calculation of the residue for evaluating the fit
#    
#    Y=Y';
#    yy=yy';
#    residue=x-Y;
#    #figure;plot(t,residue);title('Residue');xlabel('time(ps)');ylabel('Data-Fit');
#    
#    #FFt of the residue
#    
#    #Y1 = fft(residue,4096);
#    Y1=fft(residue);
#    NN=length(Y1);
#    Pyy= Y1.*conj(Y1)/(NN-1);
#    frequ = 1/(t(2)-t(1))/NN*(0:((NN/2)-1))*100/3*30;#in GHz
#    FFTResidue=Pyy(1:NN/2);
#    frequ=frequ';
#    #figure;
#    subplot(3,2,4);
#    plot(frequ,Pyy(1:NN/2));title('FFT of the Residue');xlabel('freq (GHz)');ylabel('FFT');
#    
#    #spectrum of the fit in a "Ramanlike" fashion, meaning without considering the phase or phase zero
#    
#    b=B/2/pi*100/3*30;#in GHz
#    f=W*100/3/2/pi*30;#in GHz
#    Wmax=max(W);
#    freq=0:Wmax/2/pi/1000*100/3*30:1.5*Wmax/2/pi*100/3*30;#in GHz
#    freq=freq';
#    
#    res=zeros(length(freq),length(W));
#    
#    for i=1:length(W)
#       if W(i)==0
#          res(:,i)=0;
#       else 
#          #res(:,i)=C(i)*imag(1i*exp(-1i*fi(i))*f(i)./(f(i)^2-freq.^2-2j.*freq*b(i)));
#          res(:,i)=C(i)*imag(1*f(i)./(f(i)^2-freq.^2-2j.*freq*b(i)));
#      end
#      end
#      
#      spectrum=sum(res,2);
#    #figure;
#    subplot(3,2,6);
#    plot(freq,spectrum,freq,res);title('Raman-like Data Spectrum');xlabel('freq (GHz)');ylabel('spectrum');
#    
#    
#    #for SiO2
#    t0=0.26;
#    Cp=C.*exp(B*t0);
#    fip=W*t0-fi;
#    
#    for i=1:length(W)
#          yp(:,i)=Cp(i).*exp(-B(i).*t).*cos(W(i).*t+fip(i));
#       end
#    
#    
#    
#    #Los resultados en GHz
#    tau=1./B;#in picoseconds
#    fG=W*100/3/2/pi*30;#in GHz
#    Q=W/2./B;#factor de calidad
#    
#    'Frecuencia en GHz, tau en ps, fase en grados, Amplitud, Q'
#    format shortE
#    resultsGHz=[fG 
#       tau
#       fi*180/pi
#       C
#       Q];
#    resultsGHz=resultsGHz'
#    
