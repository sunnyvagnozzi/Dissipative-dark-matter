       Implicit none
       real*8 mergedpk(1:1094),Pk(1:547),kappa(1:547),mhalo(1:547)
       real*8 dndm(1:547)
       real*8 p1,A1,f,nu1,S,acc,h,H0,rhocr,pi,GN,MP,Sk
       real*8 Omegab,Omegabp,rhob,rhobp,msun,conv,deltac,Sknew
       integer i,i1,i2,i3,i4,i5

c      Define constants

       pi = 3.14159269
       MP = 1.221E28
       GN = 1/(MP*MP)
       h = 0.7
       H0 = h*2.133E-33
       rhocr = 3.0*H0*H0/(8.0*pi*GN)
       Omegab = 0.022/(h*h)
       Omegabp = (0.14/(h*h)) - Omegab
       rhobp = Omegabp*rhocr
       rhob = Omegab*rhocr
       msun = 1.116d66
       conv = (3.08568025E29)/1.973
       p1=0.3
       A1=0.3222
       deltac=1.686
       acc=0.1

c      Open the file, which contains two columns: k and P(k)

       open(12,file="powermdmx024accurate")

c      Read the file and put it into an array

       read(12,*) (mergedpk(i),i=1,1094)

c      Split the array into even entries [P(k)] and odd entries [k]

       do 1 i1=1, 547, +1
       kappa(i1)=mergedpk(2*i1-1)
       Pk(i1)=mergedpk(2*i1)
1      continue

c      Create an array with the halo masses corresponding to k, in eV
c      Uncomment last line to get mass in units of 10^12 solar masses

       do 2 i2=1, 547, +1
       mhalo(i2)= 1.333*pi*(rhobp+rhob)*((2.5*conv/(kappa(i2)*h))**3)
c       mhalo(i2)=h*mhalo(i2)/(1.0d12*msun)
2      continue

c       Write k and M(k), this is just a check

c       do 3 i3=1, 92, +1
c       write(6,*) kappa(i3), mhalo(i3)
c3      continue

c      Compute variance

       Sk=0.0

c      Variance integral, Eq.(8) in 1412.2133

       do 4 i4=1, 547, +1
       Sk=Sk+acc*kappa(i4)*kappa(i4)*kappa(i4)*Pk(i4)/(2.0*pi*pi)
       nu1=deltac*deltac/Sk
       f=A1*(dexp(-nu1*0.5))*((2.0*nu1/pi)**(0.5))
       f=f*(1.0 + ((nu1)**(-p1)))
       dndm(i4)=(1.0/(12.0*pi*pi))*(rhobp+rhob)
       dndm(i4)=dndm(i4)*nu1*f*Pk(i4)*kappa(i4)*kappa(i4)*kappa(i4)
       dndm(i4)=dndm(i4)/(deltac*deltac*mhalo(i4))
       dndm(i4)=dndm(i4)*((conv)**(3.0))
4      continue

c      Write stellar mass in units of 10^12 solar masses and mass function

       do 5 i5=1, 547, +1
       write(6,*) dlog10(0.18*mhalo(i5)/msun), dlog10(dndm(i5)/0.18)
5      continue

       end
