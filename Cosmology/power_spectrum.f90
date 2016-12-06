           Implicit none
           real*8 MP, GN, h, H0, rhocr, conv, k, Phibegin
           real*8 f,S,f1z, f2z, f1, f2, l6, norm,rs, drs, R9, cs
           real*8 xedum,taudotdum,adum,astar1,astar2,deta,eta,etastar
           real*8 alpha2cp,mF2,me,PPip, PPipz,N3, N3z, N4, N4z, N5, N5z
           real*8 Tdec,Ip,dThrP4p,dThrP4pz,dThrP5p,dThrP5pz 
           real*8 msun,R1,dThrP3p,dThrP3pz,dThrP2p,dThrP2pz 
           real*8 dN3, dN3z, dN4, dN4z, rhogp,dN5, dN5z
           real*8 M1,dn9,ThrP4p,ThrP4pz,ThrP5p,ThrP5pz,dn1new
           real*8 dThr5p,dThr3p,dThr4p,dThr5pz,dThr3pz,dThr4pz
           real*8 n07,zeta,ThrP3p,ThrP3pz,ThrP2p,ThrP2pz 
           real*8 ThrP1p,ThrP1pz,dThrP1p,dThrP1pz
           real*8 ThrP0p,ThrP0pz,dThrP0p,nu1,deltac,dThrP0pz 
           real*8 Thr3p, Thr3pz, Thr4p, Thr4pz, Thr5p, Thr5pz
           real*8 H1b,check, dcheck,detau, dtau,tau, etap 
           real*8 Zeff,Phizbegin, a, Phi, Phiz, rhogamma, rhoneu
           real*8 rhor, rho, H1, da, dn, dPhi, dPhiz, y, y2, e
           real*8 g(0:20000),delta, ddelta, Thr0, Thr1, dThr1, dThr0
           real*8 acc2, az(0:2000),etaz(0:2000),gz(0:2000)
           real*8 deltaz, ddeltaz, Thr0z, Thr1z, dThr1z, dThr0z
           real*8 v, vz, dn6, ddeltabpz, ddeltabp, dv, dvz, T, D1new
           real*8 rhob, rhobz, N0, N0z, N1, N1z, vb, vbz, c1, c2, c3
           real*8 Omegab, deltabz, deltab, taudot,vbp,vbpz,dvbp,dvbpz
           real*8 ddeltab, ddeltabz, dvb, dvbz, R, dN1, dN1z
           real*8 s1,dN0,dN0z, Xe, Tgamma, deltaH, P,Delta2
           real*8 apeebles, bpeebles, dpeebles,vbS(0:2000),vbzS(0:2000)
           real*8 mF1,detadum,PPiS(0:2000),PPizS(0:2000)
           real*8 alpha2,alphap,alphap2,beta,dXe,c9,c4,x
           real*8 YHe, YHep,Thr0p, Thr0pz, Thr1p, Thr1pz, Tgammap
           real*8 taudotp, Rp, Xep, dThr0p, dThr0pz, dThr1p, dThr1pz
           real*8 Omegabp, alpha2c, deltabp, deltabpz,rhobp, dk
           real*8 H1a, O1, dD1, dap, ap, D1, s8,astar, n07star
           DOUBLE PRECISION J(0:7000),T1(0:7000), dC(0:7000),C(0:7000)  
           DOUBLE PRECISION z55, l99, pi, alate,Xedecoupling
           real*8 T1z(0:7000), stupid(0:7000), acc, lambda   
           real*8 ThrP1z,Thr0S(0:2000), Thr0zS(0:2000), PhiS(0:2000) 
           real*8 Thr3, Thr3z, Thr4, Thr4z, Thr5, Thr5z,PhizS(0:2000)
           real*8 dThr3, dThr3z, dThr4, dThr4z, dThr5, dThr5z
           real*8 dThrP1z, ThrP0, ThrP0z, ThrP1, p1,dThrP2, ThrP2z
           real*8 ThrP3, ThrP3z, ThrP5, ThrP5z, ThrP4, ThrP4z
           real*8 PPiz, PPi,dThrP0, dThrP0z, dThrP1, A1,dThrP2z
           real*8 dThrP3, dThrP3z, dThrP5, dThrP5z, dThrP4, dThrP4z
           real*8 Thr2, Thr2z, Psi, Psiz, PsiS(0:2000), PsizS(0:2000)
           real*8 dThr2, dThr2z, dPsi, dPsiz, N2, N2z, dN2, dN2z   
           real*8 dOM, dI1, dI2, I1, I2, dNeff, t8, ThrP2,J88
           real*8 Neff,Thr2p, Thr2pz,dThr2p, dn7,dThr2pz
           Integer*4 k4,k7,l,j9,i,l4,m,npeak2,npeak,n,k8,m6,N6

c all equations refenced are from Dodelson's book, except when stated otherwise

           open(unit=4,file="fig1resxacc")
    
c set some constants

           Neff = 3.046 ! neglect energy density of exotic states, not very important for power spectrum. 

           alate = 2.0d-2  ! alate is same notation as dodelson..(see around eq.7.4)
           pi = 3.14159269   
           me = 0.511d6 ! electron mass in eV
           MP = 1.221E28 ! Mplank mass in eV  
           GN = 1/(MP*MP)
           conv = (3.08568025E29)/1.973  ! converts Mpc to 1/eV.
           msun = 1.116d66


           h = 0.7 ! h is usual hubble coefficient
           H0 = h*2.133E-33 ! Hubble constant units of eV, eq 1.5.
           rhocr = 3.0*H0*H0/(8.0*pi*GN) ! critical density units of eV^4, Eq.1.3

           Omegab = 0.022/(h*h) ! ordinary baryon part.
           Omegabp = (0.14/(h*h)) - Omegab  !  cold dark matter part

c           write(6,*) 'Omegab=',Omegab, 'Omegabp =', Omegabp

c set Y helium mass fraction

           YHe = 0.24

c set mirror helium mass fraction

           YHep = 0 ! for dark sector model, only two components, F1 and F2

c set x value...
 
           x = 0.24 ! T'/T

c set hidden sector parameters

          alphap = 1.0/137.0
          alphap2 = alphap*alphap
          mF1 = me
          mF2 = 940.0d6  ! set equal to proton mass for now
          Zeff = 1.0
          alpha2c = 9.78/(137.0*137.0*me*me) ! coefficient of Eq.12 in arXiv:1208.6022
          alpha2cp = 9.78*alphap*alphap/(mF1*mF1) ! same as alpha2c above but for dark sector


c          work out conformal time today

           rhocr = 3.0*H0*H0/(8.0*pi*GN) ! units of eV^4, Eq.1.3

           ap = 1.0d-8
           n07 = 0
           D1 = 0
6          O1 = Omegab + Omegabp
           rhob = Omegab*rhocr/(ap*ap*ap)
           rhobp = Omegabp*rhocr/(ap*ap*ap)
           rhogamma = rhocr*2.47E-5/(h*h*ap*ap*ap*ap) ! Eq.2.70
           rhogp = rhogamma*x*x*x*x
           rhoneu = rhogamma*(Neff*7.0/8.0)*((4.0/11.0)**(1.33333)) !  Eq.2.77
           rhor = rhogp+rhogamma+rhoneu

           rho = rhor + rhob + rhobp
           rho = rho + (1.0-O1)*rhocr

           H1a = (8.0*pi*GN*0.333333*rho)**(0.5)

           dn9 = 0.0001*ap/(H1a*ap*ap)                   ! Eq.2.41
           n07 = n07 + dn9
           if(ap.gt.alate)then
           dD1 = 2.5*O1*0.0001*ap/((ap*H1a/H0)**(3.0)) ! Eq.7.77
           D1 = D1 + dD1
           endif
           ap = ap*1.0001
           if(ap.lt.1.0)then
           goto 6
           endif

           D1=D1*H1a/H0

           write(6,*) 'n0 (Mpc) =',n07/conv,'D1=',D1

c set initial k value

           S = 0 ! Eq.8 of Schneider paper.
           
           k = 0.001  ! starting value of k units h Mpc^-1  
           
           k = h*k/conv  ! converts k in (h Mpc^-1 to eV), so from now on k has units eV
           
c start the k loop ...work out power spectrum for k = 0.001 till k = 1.0
           
3          a = 3.0d-6  ! initial value of a used in solving the equations 
           if((k*conv/h).lt.0.1)then
           a = 8.0d-6  ! can have large a initially for smaller k...see figure 7.3 of dodelson 
           endif

c as we make this initial value smaller..our results for power spectrum shouldn't change....
           
c       initial condition
c first for ordinary particles
        
           Xep = 1.0
           Xe = 1.0
        
           Phibegin = (1.0E-12)/(k**1.5)  ! arbitrary initial value, value not important as rescale to get correct power spectrum on large scales 

           Phizbegin = 0           
     
           Phi = Phibegin
           Phiz = Phizbegin   ! Phiz is imaginary part, Phi is real part..same for all other variables   
           
           Thr0 = 0.5*Phibegin  ! Eq.6.12
           Thr0z = 0.5*Phizbegin
           Thr0p = 0.5*Phibegin
           Thr0pz = 0.5*Phizbegin
           
           N0 = Thr0 ! Eq.6.11
           N0z = Thr0z
           
           deltab = 1.5*Phibegin  ! Eq.6.15,6.12
           deltabz = 1.5*Phizbegin

c need H for initial condition, Eq.6.16

           rhob = Omegab*rhocr/(a*a*a)
           rhobp = Omegabp*rhocr/(a*a*a)
           rhogamma = rhocr*2.47E-5/(h*h*a*a*a*a) ! Eq.2.70
           rhogp = rhogamma*x*x*x*x
           rhoneu = rhogamma*(Neff*7.0/8.0)*((4.0/11.0)**(1.33333)) !  Eq.2.77
           rhor = rhogp+rhogamma+rhoneu
           rho = rhor + rhob + rhobp
           H1 = (8.0*pi*GN*0.333333*rho)**(0.5)

           Thr1 = -k*Phibegin/(6.0*a*H1)     ! Eq.6.16
           Thr1z = -k*Phizbegin/(6.0*a*H1)

           Thr2 = 0
           Thr2z = 0

           N1 = Thr1  ! Eq.6.16
           N1z = Thr1z

           vb = -k*Phizbegin/(2.0*a*H1) ! Eq.6.16
           vbz = k*Phibegin/(2.0*a*H1)


c now for initial condition for mirror particles, denoted with supscript p (for prime)

           deltabp = deltab
           deltabpz = deltabz
           
           Thr1p = -k*Phibegin/(6.0*a*H1)     ! Eq.6.16
           Thr1pz = -k*Phizbegin/(6.0*a*H1)
             
              
           vbp = vb
           vbpz = vbz
          
           Psi = -Phi
           Psiz = -Phiz
           
           Thr2p = 0
           Thr2pz = 0

c set intial conditons for all higher moments to zero
          
           N2 = 0
           N2z = 0
           N3 = 0
           N3z = 0
           N4 = 0
           N4z = 0
           N5 = 0
           N5z = 0

c add it extra l to l=5
 
           Thr3 = 0 
           Thr3z = 0 
           Thr3p = 0 
           Thr3pz = 0 
           Thr4 = 0 
           Thr4z = 0 
           Thr4p = 0 
           Thr4pz = 0 
           Thr5 = 0 
           Thr5z = 0 
           Thr5p = 0 
           Thr5pz = 0 
           ThrP0 = 0 
           ThrP0z = 0 
           ThrP0p = 0 
           ThrP0pz = 0 
           ThrP1 = 0 
           ThrP1z = 0 
           ThrP1p = 0 
           ThrP1pz = 0 
           ThrP2 = 0 
           ThrP2z = 0 
           ThrP2p = 0 
           ThrP2pz = 0 
           ThrP3 = 0 
           ThrP3z = 0 
           ThrP3p = 0 
           ThrP3pz = 0 
           ThrP4 = 0 
           ThrP4z = 0 
           ThrP4p = 0 
           ThrP4pz = 0 
           ThrP5 = 0 
           ThrP5z = 0 
           ThrP5p = 0 
           ThrP5pz = 0 
            
           PPi = 0
           PPiz = 0
           PPip = 0
           PPipz = 0

c finished initial conditions

c         start time loop              
                
        
1          rhob = Omegab*rhocr/(a*a*a)
           rhobp = Omegabp*rhocr/(a*a*a)
           
           rhogamma = rhocr*2.47E-5/(h*h*a*a*a*a)
           rhogp = rhogamma*x*x*x*x 
           rhoneu = rhogamma*(Neff*7.0/8.0)*((4.0/11.0)**(1.33333))
           rhor = rhogp+rhogamma+rhoneu

           rho = rhor + rhob + rhobp ! rho is total matter + radiation density
           
           H1 = (8.0*pi*GN*0.333333*rho)**(0.5)    ! Eq.1.3  for times of interest can negect cosmological constant


           dn = min(a*a*1.2E33,0.001*alate*alate*1.0E33) 
           if((k*conv/h).gt.1.8)then
           dn = dn*0.1
           endif

           da = H1*a*a*dn  ! Eq.2.41  ...n is neta..conformal time.
           
c         Phi equation 5.27  (z is imaginary part)

         
           s1 = 4.0*pi*GN*a*a
           
           dPhi=s1*(rhobp*deltabp+4.0*rhogamma*Thr0)
           dPhi=dPhi+s1*(rhob*deltab+4.0*rhoneu*N0+4.0*rhogp*Thr0p) !  Thr0 = N0 if neglect  baryons       
           dPhi = dPhi - k*k*Phi
           dPhi = (dPhi*0.33333/(a*H1))*dn + a*H1*Psi*dn

           dPhiz=s1*(rhob*deltabz+rhobp*deltabpz)
           dPhiz=dPhiz+s1*4.0*(rhogamma*Thr0z+rhoneu*N0z+rhogp*Thr0pz)
           dPhiz = dPhiz- k*k*Phiz
           dPhiz = (dPhiz*0.33333/(a*H1))*dn + a*H1*Psiz*dn
           
c         difusion equation for Psi

           Psi =-Phi - (8.0*s1*(rhogamma*Thr2+rhoneu*N2)/(k*k)) !  Eq.5.33
           Psi = Psi - (8.0*s1*rhogp*Thr2p/(k*k))           
           Psiz =-Phiz-(8.0*s1*(rhogamma*Thr2z+rhoneu*N2z)/(k*k))
           Psiz = Psiz - (8.0*s1*rhogp*Thr2pz/(k*k))           
                        
c         Thetar equation 7.11
                       
           dThr0 = -k*Thr1*dn - dPhi
           
           dThr0z = - k*Thr1z*dn - dPhiz
           
           dN0 = -k*N1*dn - dPhi
           
           dN0z = - k*N1z*dn - dPhiz

           dThr0p = -k*Thr1p*dn - dPhi
           
           dThr0pz = - k*Thr1pz*dn - dPhiz
                   
           
c          eq.712

ccccccccccccccc mirror dark matter equations
cccccccccccccccccccccccccccccccccccccccccccc
           if(x.lt.0.01)then  
           taudotp = 0  ! same as cold dark matter
           goto 45
           endif

           Tgammap = (15.0*rhogp/(pi*pi))**(0.25) ! temperature of mirror particles
           Ip = alphap2*Zeff*Zeff*mF1*0.5   ! ionization energy of F1 bound state
           zeta = 40.0 ! might need to work it out more accurately, use 40 for now Eq.44 of 1409.7174
           Tdec = Ip/zeta ! this is approx temperature' of decoupling

           if(Tgammap.gt.(1.5*Tdec))then  ! if T' > 2.0Tdec means don't need to solve equations, mirror states are ionized  
           goto 45
           endif
           if(Tgammap.lt.(0.20*Tdec))then ! if T < 0.20Tdec, don't need to solve equations, mirror states fully combined.
           goto 45
           endif

c need smaller time step to solve these rate equations. Use dn7 = dn/100.0 and do 100 times
           
           dn7 = dn/100.0
           do j9=1,100,+1
           alpha2 = alpha2cp*(Ip/Tgammap)**(0.5)
           alpha2 = alpha2*log(Ip/Tgammap)
           beta = alpha2*((mF1)**1.5)
           beta = beta * ((Tgammap/6.283)**1.5)
           beta = beta*dexp(-Ip/Tgammap)
           c4 = (1.0-YHep)*(rhobp/mF2)*alpha2
           
           dXe = (1.0-Xep)*beta - Xep*Xep*c4   ! equations here and below given in Eq.12 of 1208.6022
           if(Xep.lt.0.99)then
           apeebles = ((3.0*Ip)**3)*H1/(64.0*pi*pi)
           apeebles=apeebles/((1.0-Xep)*(1.0-YHep)*(rhobp/mF2))
           bpeebles = 8.227*6.582E-16 ! need to scale this...not sure how two photon decay scales with alpha' and mF1...check peebles 1968 paper.
           dpeebles = beta*exp(0.75*Ip/Tgammap)
           c9 = (apeebles+bpeebles)/(apeebles + bpeebles + dpeebles)
           else
           c9 = 1.0d-2
           endif
           dXe = dXe*c9
           dXe = dXe*a*dn7
           Xep = Xep+dXe         
           enddo
           
45         taudotp = -(0.0692/(a*a))*Xep*(1.0-YHep)*Omegabp*h*H0

           if(Tgammap.gt.0.5)then
           taudotp = -(0.0692/(a*a))*Xep*(1.0-0.5*YHep)*Omegabp*h*H0 ! effects of helium
           endif


           Rp = 0.75*rhobp/rhogp

           dvbp=-a*H1*vbp*dn+k*Psiz*dn+dn*(taudotp/Rp)*(vbp-3.0*Thr1pz)
           
           dvbpz=-a*H1*vbpz*dn-k*Psi*dn+dn*(taudotp/Rp)*(vbpz+3.0*Thr1p)

           dThr1p=0.3333*k*dn*(Thr0p+Psi-2.0*Thr2p)
           dThr1p = dThr1p+dn*taudotp*(Thr1p+vbpz*0.33333)

           dThr1pz=0.3333*k*dn*(Thr0pz+Psiz-2.0*Thr2pz)
           dThr1pz = dThr1pz +dn*taudotp*(Thr1pz-vbp*0.33333)

           dThr2p=k*dn*(0.4*Thr1p-0.6*Thr3p)
           dThr2p=dThr2p+taudotp*dn*(Thr2p-(PPip*0.1))

           dThr2pz=k*dn*(0.4*Thr1pz - 0.6*Thr3pz)
           dThr2pz=dThr2pz+taudotp*dn*(Thr2pz-(PPipz*0.1))

           dThr3p = k*dn*(0.4286*Thr2p - 0.5714*Thr4p)+taudotp*dn*Thr3p
           dThr3pz=k*dn*(0.4286*Thr2pz- 0.5714*Thr4pz)+taudotp*dn*Thr3pz

           dThr4p = k*dn*(0.4444*Thr3p - 0.5555*Thr5p)+taudotp*dn*Thr4p
           dThr4pz =k*dn*(0.4444*Thr3pz-0.5555*Thr5pz)+taudotp*dn*Thr4pz

           dThr5p = k*dn*(0.4545*Thr4p)+taudotp*dn*Thr5p
           dThr5pz = k*dn*(0.4545*Thr4pz)+taudotp*dn*Thr5pz

           dThrP0p = -k*ThrP1p*dn + taudotp*dn*(ThrP0p - PPip*0.5)
           dThrP0pz = -k*ThrP1pz*dn + taudotp*dn*(ThrP0pz - PPipz*0.5)

           dThrP1p = k*dn*(0.3333*ThrP0p - 0.6667*ThrP2p) 
           dThrP1p = dThrP1p +taudotp*dn*ThrP1p 

           dThrP1pz=k*dn*(0.3333*ThrP0pz- 0.6667*ThrP2pz)
           dThrP1pz = dThrP1pz+taudotp*dn*ThrP1pz 

           dThrP2p = k*dn*(0.4*ThrP1p - 0.6*ThrP3p) 
           dThrP2p = dThrP2p +taudotp*dn*(ThrP2p-PPip*0.1) 

           dThrP2pz=k*dn*(0.4*ThrP1pz-0.6*ThrP3pz)
           dThrP2pz = dThr2pz +taudotp*dn*(ThrP2pz-PPipz*0.1) 

           dThrP3p=k*dn*(0.4286*ThrP2p-0.5714*ThrP4p)+taudotp*dn*ThrP3p 
           dThrP3pz =k*dn*(0.4286*ThrP2pz-0.5714*ThrP4pz)
           dThrP3pz=dThrP3pz+taudotp*dn*ThrP3pz 

           dThrP4p=k*dn*(0.4444*ThrP3p-0.5555*ThrP5p)+taudotp*dn*ThrP4p 
           dThrP4pz =k*dn*(0.4444*ThrP3pz-0.5555*ThrP5pz)
           dThrP4pz = dThrP4pz+taudotp*dn*ThrP4pz 

           dThrP5p = k*dn*0.4545*ThrP4p+taudotp*dn*ThrP5p 
           dThrP5pz = k*dn*0.4545*ThrP4pz+taudotp*dn*ThrP5pz 
           

 
           ddeltabp = + k*vbpz*dn - 3.0*dPhi
           
           ddeltabpz = -k*vbp*dn - 3.0*dPhiz 
            
          
c  here
           
cccccccccccc   ordinary baryon equations cccccccccccccccccccccccccccc           

           PPi = Thr2 + ThrP2 + ThrP0
           PPiz = Thr2z + ThrP2z + ThrP0z
           PPip = Thr2p + ThrP2p + ThrP0p
           PPipz = Thr2pz + ThrP2pz + ThrP0pz
         
           
           Tgamma = (15.0*rhogamma/(pi*pi))**(0.25)
           if(Tgamma.gt.0.40)then ! don't need to solve rate equations for T > 0.4 eV...as H fully ionized
           goto 46
           endif
           if(Tgamma.lt.0.08)then ! don't need to solve rate equations for T < 0.08 eV as H fully combined
           goto 46
           endif
           
c need smaller time step to solve these rate equations. Use dn7 = dn/100.0 and do 100 times
           
           dn7 = dn/100.0
           do j9=1,100,+1
           alpha2 = alpha2c*((13.6/Tgamma)**(0.5))
           alpha2 = alpha2*log(13.6/Tgamma)
           beta = alpha2*((0.5*Tgamma*0.511E6/pi)**1.5)
           beta = beta*(dexp(-13.6/Tgamma))
           c4 = (1.0-YHe)*(rhob/9.4E8)*alpha2
           

   
           dXe = (1.0-Xe)*beta - Xe*Xe*c4
           if(Xe.lt.0.99)then
           apeebles = ((3.0*13.6)**3)*H1
           apeebles= apeebles/(64.0*pi*pi*(1.0-Xe)*(1-YHe)*(rhob/9.4E8))
           bpeebles = 8.227*6.582E-16
           dpeebles = beta*exp(0.75*13.6/Tgamma)
           c9 = (apeebles+bpeebles)/(apeebles + bpeebles + dpeebles)
           else
           c9 = 1.0d-2
           endif
           dXe = dXe*c9
           dXe = dXe*a*dn7
                
           Xe = Xe+dXe         
           enddo
           
46         taudot = -(0.0692/(a*a))*Xe*(1.0-Yhe)*Omegab*h*H0  ! Eq.4.62 and Eq.3.45, assumes Helium is combined

           if(Tgamma.gt.0.5)then
           taudot = -(0.0692/(a*a))*Xe*(1.0-0.5*Yhe)*Omegab*h*H0  ! assumes Helium is singly ionized so have slightly more electrons.
           endif
           
           
           R = 0.75*rhob/rhogamma
         
           dThr1=0.33333*k*dn*(Thr0+Psi-2.0*Thr2)
           dThr1=dThr1+dn*taudot*(Thr1+vbz*0.33333)

           dThr1z=0.33333*k*dn*(Thr0z+Psiz-2.0*Thr2z)
           dThr1z = dThr1z +dn*taudot*(Thr1z-vb*0.33333)
           
c           put in difusion...quadrupole term           

           
           
           dThr2 = k*dn*(0.4*Thr1 - 0.6*Thr3)
           dThr2=dThr2+taudot*dn*(Thr2 - (PPi*0.1))
           dThr2z=k*dn*(0.4*Thr1z - 0.6*Thr3z)
           dThr2z=dThr2z+taudot*dn*(Thr2z-(PPiz*0.1))
           dThr3 = k*dn*(0.4286*Thr2 - 0.5714*Thr4)+taudot*dn*Thr3
           dThr3z = k*dn*(0.4286*Thr2z- 0.5714*Thr4z)+taudot*dn*Thr3z
           dThr4 = k*dn*(0.4444*Thr3 - 0.5555*Thr5)+taudot*dn*Thr4
           dThr4z = k*dn*(0.4444*Thr3z - 0.5555*Thr5z)+taudot*dn*Thr4z
           dThr5 = k*dn*(0.4545*Thr4)+taudot*dn*Thr5
           dThr5z = k*dn*(0.4545*Thr4z)+taudot*dn*Thr5z
           dThrP0 = -k*ThrP1*dn + taudot*dn*(ThrP0 - PPi*0.5)
           dThrP0z = -k*ThrP1z*dn + taudot*dn*(ThrP0z - PPiz*0.5)
           dThrP1 = k*dn*(0.3333*ThrP0 - 0.6667*ThrP2) +taudot*dn*ThrP1 
           dThrP1z=k*dn*(0.3333*ThrP0z- 0.6667*ThrP2z)+taudot*dn*ThrP1z 
           dThrP2 = k*dn*(0.4*ThrP1 - 0.6*ThrP3)
           dThrP2=dThrP2+taudot*dn*(ThrP2-PPi*0.1)
           dThrP2z=k*dn*(0.4*ThrP1z-0.6*ThrP3z)
           dThrP2z=dThrP2z+taudot*dn*(ThrP2z-PPiz*0.1) 
           dThrP3 = k*dn*(0.4286*ThrP2-0.5714*ThrP4)+taudot*dn*ThrP3 
           dThrP3z =k*dn*(0.4286*ThrP2z-0.5714*ThrP4z)+taudot*dn*ThrP3z 
           dThrP4 = k*dn*(0.4444*ThrP3-0.5555*ThrP5)+taudot*dn*ThrP4 
           dThrP4z =k*dn*(0.4444*ThrP3z-0.5555*ThrP5z)+taudot*dn*ThrP4z 
           dThrP5 = k*dn*0.4545*ThrP4+taudot*dn*ThrP5 
           dThrP5z = k*dn*0.4545*ThrP4z+taudot*dn*ThrP5z 
             
           dN2 = k*dn*(0.4*N1 - 0.6*N3)
           dN2z = k*dn*(0.4*N1z - 0.6*N3z)
           dN3 = k*dn*(0.4286*N2 - 0.5714*N4)
           dN3z = k*dn*(0.4286*N2z - 0.5714*N4z)
           dN4 = k*dn*(0.4444*N3 - 0.5555*N5)
           dN4z = k*dn*(0.4444*N3z - 0.5555*N5z)
           dN5 = k*dn*(0.4545*N4)
           dN5z = k*dn*(0.4545*N4z)
           
           
           
c          put in difusion            
           
           
          
           dN1 = 0.33333*k*dn*(N0+Psi-2.0*N2) 
           
           dN1z = 0.33333*k*dn*(N0z + Psiz - 2.0*N2z)
           
           
           
c          eq.7.13

           
           ddeltab = + k*vbz*dn - 3.0*dPhi
           
           ddeltabz = -k*vb*dn - 3.0*dPhiz
           
    
           
c          eq. 7.14

c dv = -a*H1*v*dn + k*Psiz*dn
c dvz = -a*H1*vz*dn - k*Psi*dn
           
           
           dvb = -a*H1*vb*dn + k*Psiz*dn + dn*(taudot/R)*(vb -3.0*Thr1z)
           
           dvbz = -a*H1*vbz*dn - k*Psi*dn +dn*(taudot/R)*(vbz+3.0*Thr1)
           
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc          
                                 
           
           Phi = Phi + dPhi
           Phiz = Phiz + dPhiz
           
           Psi = Psi + dPsi
           Psiz = Psiz + dPsiz
           
           
           N0 = N0 + dN0
           N0z = N0z + dN0z
           
           N1 = N1 + dN1
           N1z = N1z +dN1z
           
           N2 = N2 + dN2
           N2z = N2z +dN2z
           
           N3 = N3 + dN3
           N3z = N3z +dN3z
           N4 = N4 + dN4
           N4z = N4z +dN4z
           N5 = N5 + dN5
           N5z = N5z +dN5z
           
           Thr0 = Thr0 + dThr0
           Thr0z = Thr0z + dThr0z
           Thr0p = Thr0p + dThr0p
           Thr0pz = Thr0pz + dThr0pz
           
           Thr1 = Thr1 + dThr1
           Thr1z = Thr1z + dThr1z
           Thr1p = Thr1p + dThr1p
           Thr1pz = Thr1pz + dThr1pz
          
           Thr2 = Thr2 + dThr2
           Thr2z = Thr2z + dThr2z  
           Thr2p = Thr2p + dThr2p
           Thr2pz = Thr2pz + dThr2pz  

c extra stuff here
           
           Thr3 = Thr3 + dThr3
           Thr3z = Thr3z + dThr3z
           Thr3p = Thr3p + dThr3p
           Thr3pz = Thr3pz + dThr3pz
           Thr4 = Thr4 + dThr4
           Thr4z = Thr4z + dThr4z
           Thr4p = Thr4p + dThr4p
           Thr4pz = Thr4pz + dThr4pz
           Thr5 = Thr5 + dThr5
           Thr5z = Thr5z + dThr5z
           Thr5p = Thr5p + dThr5p
           Thr5pz = Thr5pz + dThr5pz
           ThrP0 = ThrP0 + dThrP0
           ThrP0z = ThrP0z + dThrP0z
           ThrP0p = ThrP0p + dThrP0p
           ThrP0pz = ThrP0pz + dThrP0pz
           ThrP1 = ThrP1 + dThrP1
           ThrP1z = ThrP1z + dThrP1z
           ThrP1p = ThrP1p + dThrP1p
           ThrP1pz = ThrP1pz + dThrP1pz
           ThrP2 = ThrP2 + dThrP2
           ThrP2z = ThrP2z + dThrP2z
           ThrP2p = ThrP2p + dThrP2p
           ThrP2pz = ThrP2pz + dThrP2pz
           ThrP3 = ThrP3 + dThrP3
           ThrP3z = ThrP3z + dThrP3z
           ThrP3p = ThrP3p + dThrP3p
           ThrP3pz = ThrP3pz + dThrP3pz
           ThrP4 = ThrP4 + dThrP4
           ThrP4z = ThrP4z + dThrP4z
           ThrP4p = ThrP4p + dThrP4p
           ThrP4pz = ThrP4pz + dThrP4pz
           ThrP5 = ThrP5 + dThrP5
           ThrP5z = ThrP5z + dThrP5z
           ThrP5p = ThrP5p + dThrP5p
           ThrP5pz = ThrP5pz + dThrP5pz

c  end extra stuff
           

           deltab = deltab + ddeltab
           deltabz = deltabz + ddeltabz
            
           deltabp = deltabp + ddeltabp
           deltabpz = deltabpz + ddeltabpz
           
           
           vb = vb +dvb
           vbz = vbz + dvbz
           
           vbp = vbp +dvbp
           vbpz = vbpz + dvbpz
           
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc           
 
c work out transfer function    
            
               
           a = a + da

           if(a.lt.alate)then
           goto 1
           endif    

c finished time loop 

           
c work out transfer function    
                             
           T = (Phi*Phi + Phiz*Phiz)**(0.5)
           T = T/(Phibegin*Phibegin+Phizbegin*Phizbegin)**(0.5)    
           T = T/0.90
                       
           
           deltaH = (4.6E-5) * (0.30/(Omegab+Omegabp))
           
           P = 2.0*pi*pi*deltaH*deltaH*k*T*T/(H0**4)
c           P = P*((D1/0.76947107)**(2.0))

           Delta2 = k*k*k*P/(2.0*pi*pi)
           
           acc = 0.1
           if((k*conv/h).gt.0.06)then
           acc = 0.1
           endif

c work out S(R)

           S = S+ k*k*acc*k*P/(2.0*pi*pi) ! acc*k = dk

           deltac = 1.686

           D1new=D1/0.76947107

           nu1 = deltac*deltac/(S*D1new*D1new)
c           write(6,*)'S=',S

           p1 = 0.3
           A1 = 0.3222
           f = A1*(dexp(-nu1*0.5))*((2.0*nu1/pi)**(0.5))
           f = f*(1.0  + ((nu1)**(-p1)))

           rhobp = Omegabp*rhocr   ! need to divide by (a*a*a)
           rhob = Omegab*rhocr   ! need to divide by (a*a*a)
           dn1 = (1.0/(12.0*pi*pi))*(rhobp+rhob)
           dn1 = dn1*nu1*f*P*k*k*k/(deltac*deltac)
           M1 = 1.333*pi*(rhobp+rhob)*((2.5/k)**3) ! used c=2.5,k has units ev, rho has units ev^4, therefore M1 has units eV
           dn1 = dn1/M1  ! dn has units eV^3
           dn1 = dn1*D1new*D1new
           dn1 = dn1*((conv/h)**3) ! converts dn to units (h/Mpc)^3

           P = P*(h/conv)**3 ! convert P into Mpc/h ^3
           write(6,*) k*conv/h, P
c           write(6,*) h*M1/(1.0d12 *msun), dn1
c           write(4,*) k*conv/h, h*M1/msun,P, dn1
c           write(6,*) k*conv/h, dn1
c           write(6,*) (rhobp+rhob)*(conv**(3.0))/(1.0d10*msun)


           k = k + acc*k
           
           if((k*conv/h).lt.6.0)then
           goto 3
           endif

           end
