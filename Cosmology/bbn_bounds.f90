       Implicit none
       real*8 Tg, Tgp, Tn, dTg, dTgp, dTn, Tinitial, Tabc
       real*8 Rg, Pg, Re, Pe, Rgp, Pgp, Rep, Pep, Rn
       real*8 dRe, Dpe, dRep, dPep
       real*8 ep, pi, me, mf, mp, alpha, z, s, ds, t, dt, tp
       real*8 x, xp, d, c, u, up, du
       real*8 ne, dne, E, dE, Eplus, dEplus
       real*8 sig, sigvmol1, dsigvmol1, sigvmol, dsigvmol, rhs
       real*8 K2, b, a, dK2, da, depsilon, Yp, Ype0
       real*8 f1, df1, f1p, df1p, f2, df2, f2p, df2p
       real*8 r, n, o, p, np, op, pp
       real*8 Xn, a1, a2, a3, a4, dXn, da1, da2, da3, da4
       real*8 lam1, lam2, lam3, lam4, dlam1, dlam2, dlam3, dlam4
       real*8 En, Ee, Awk, Q, Gwk, ga, cab, ext3, ext4
       integer i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i, j
 
       ep=4.71d-9
       do 93 j=1,1,+1
       depsilon=0
       pi=3.1415927
       me=0.511 !MeV
       mf=0.3 !MeV
       mp=1.22d22 !MeV
       alpha=1/(137.0)
       Q=1.293 !MeV
       Gwk=1.16d-11 !MeVE-2
       ga=1.257
       cab=0.9745
       ext3=((Q*Q-me*me)**(0.5))
       ext4=Q+me
       Awk=(Gwk*Gwk*(1+3*ga*ga)*cab*cab)/(2*pi*pi*pi)
       
c set initial conditions and start time loop

       Tinitial=1000.0
       Tg=Tinitial
       Tn=Tinitial
       Tgp=10.0
       Xn=0.46
       i=0
       dt=(1.5E12)/(184000000.0)
       t=0
       do 6 i6=1,1975000,+1

c choose z, maximum between me and mf
       
       if(mf.gt.me) then
       z=mf
       else
       z=me
       endif

c start calculation of RHS of 2nd law thermodynamics
c start calculation of ne

       x=me/Tg
       xp=mf/Tgp
       rhs=0

       ne=0
       dE=Tg*14/1500.0
       E=dE+me
       do 1 i1=1,1500,+1
       dne=dE*E*((E*E-me*me)**(0.5))
       dne=dne/(pi*pi*(1+exp(E/Tg)))
       ne=ne+dne
       E=E+dE
1      continue

c end calculation of ne
c start calculation of sigvmol, split in 2 integrals
c start calculation of sigvmol1

       sigvmol=0
       s=4.0*z*z
       do 3 i3=1,300,+1
       Eplus=(s)**(0.5)
       dEplus=Tg*14.0/300.0
       sigvmol1=0
       do 2 i2=1,300,+1
       c=Eplus*Eplus
       c=c/s
       c=c-0.999
       dsigvmol1=dEplus*Eplus
       dsigvmol1=dsigvmol1*((c)**(0.5))
       dsigvmol1=dsigvmol1*(exp(-Eplus/Tg))
       sigvmol1=sigvmol1+dsigvmol1
       Eplus=Eplus+dEplus
2      continue
       
c end calculation of sigvmol1
c continue calculation of sivgmol

       ds=Tg*Tg*99.0/300.0
       sig=4.0*pi*alpha*alpha*ep*ep*0.333333
       d=s-3.999999*mf*mf
       d=(d/(s-3.999999*me*me))**(0.5)
       sig=sig*(s*s+2.0*(me*me+mf*mf)*s+4.0*me*me*mf*mf)
       sig=sig*d
       sig=sig/(s*s*s)
       dsigvmol=ds*sig*sigvmol1
       dsigvmol=dsigvmol*(s-4.0*me*me)
       dsigvmol=dsigvmol*((s)**(0.5))
       sigvmol=sigvmol+dsigvmol
       s=s+ds
3      continue

c start calculation of K2

       b=me/Tg
       if(Tg.gt.100.0) then
       K2=2.0/(b*b)
       else
       a=1.0
       da=0.01
       K2=0
       do 4 i4=1,290000,+1
       b=me/Tg
       dK2=b*b*0.33333*da
       dK2=dK2*((a*a-0.999)**(1.5))
       dK2=dK2*(exp(-b*a))
       K2=K2+dK2
       a=a+da
4      continue
       endif

c continue calculation of sigvmol

       sigvmol=(0.8*sigvmol)/(8.0*me*me*me*me*Tg*Tg*K2*K2)

c end calculation of sigvmol
c continue calculation of RHS of 2nd law of thermodynamics

       rhs=ne*ne*sigvmol
       rhs=rhs/(Tn*Tn*Tn)

c end calculation of RHS of 2nd law of thermodynamics
c calculate Rg, Pg, Rgp, Pgp, Rn

       Rg=(pi*pi*Tg*Tg*Tg*Tg)/15.0
       Pg=(pi*pi*Tg*Tg*Tg*Tg)/45.0
       Rgp=(pi*pi*Tgp*Tgp*Tgp*Tgp)/15.0
       Pgp=(pi*pi*Tgp*Tgp*Tgp*Tgp)/45.0
       Rn=(7.0*pi*pi*Tn*Tn*Tn*Tn)/40.0

c calculate Re, Pe, Rep, Pep, f1, f1p, f2, f2p

       Re=0
       Pe=0
       Rep=0
       Pep=0
       f1=0
       f1p=0
       f2=0
       f2p=0
       du=0.003
       u=x+du
       up=xp+du
       do 5 i5=1,5000,+1
       dRe=u*u*du
       dRe=dRe*((u*u-x*x)**(0.5))
       dRe=dRe/(1.0+dexp(u))
       dRe=(2.0*Tg*Tg*Tg*Tg*dRe)/(pi*pi)
       dRep=up*up*du
       dRep=dRep*((up*up-xp*xp)**(0.5))
       dRep=dRep/(1+dexp(up))
       dRep=(2.0*Tgp*Tgp*Tgp*Tgp*dRep)/(pi*pi)
       dPe=du*((u*u-x*x)**(1.5))
       dPe=dPe/(1.0+dexp(u))
       dPe=(2.0*Tg*Tg*Tg*Tg*dPe)/(3.0*pi*pi)
       dPep=du*((up*up-xp*xp)**(1.5))
       dPep=dPep/(1.0+dexp(up))
       dPep=(2.0*Tgp*Tgp*Tgp*Tgp*dPep)/(3.0*pi*pi)
       df1=u*u*du
       df1=df1*((u*u-x*x)**(-0.5))
       df1=df1/(1.0+dexp(u))
       df1=(2.0*me*me*df1)/(pi*pi*Tg*Tg)
       df1p=up*up*du
       df1p=df1p*((up*up-xp*xp)**(-0.5))
       df1p=df1p/(1.0+dexp(up))
       df1p=(2.0*mf*mf*df1p)/(pi*pi*Tgp*Tgp)
       df2=du*((u*u-x*x)**(0.5))
       df2=df2/(1.0+dexp(u))
       df2=(2.0*me*me*df2)/(pi*pi*Tg*Tg)
       df2p=du*((up*up-xp*xp)**(0.5))
       df2p=df2p/(1.0+dexp(up))
       df2p=(2.0*mf*mf*df2p)/(pi*pi*Tgp*Tgp)
       Re=Re+dRe
       Pe=Pe+dPe
       Rep=Rep+dRep
       Pep=Pep+dPep
       f1=f1+df1
       f1p=f1p+df1p
       f2=f2+df2
       f2p=f2p+df2p
       u=u+du
       up=up+du
5      continue

c linearize and change Tn, Tg, Tgp

       r=8.0*pi*0.3333
       r=r/(mp*mp)
       r=r*(Rg+Re+Rn+Rgp+Rep)
       dTn=-Tn*dt*((r)**(0.5))
       n=(((3.0*(Rg+Pg+Re+Pe))/(Tg*Tg*Tg*Tg)))+f1+f2
       o=(3.0*(Rg+Pg+Re+Pe))/(Tg*Tg*Tg*Tn)
       p=rhs*((Tn/Tg)**(3.0))
       np=(((3.0*(Rgp+Pgp+Rep+Pep))/(Tgp*Tgp*Tgp*Tgp)))+f1p+f2p
       op=(3.0*(Rgp+Pgp+Rep+Pep))/(Tgp*Tgp*Tgp*Tn)
       pp=rhs*((Tn/Tgp)**(3.0))
       dTg=(o*dTn-p*dt)
       dTg=dTg/n
       dTgp=(op*dTn+pp*dt)
       dTgp=dTgp/np
       Tg=Tg+dTg
       Tgp=Tgp+dTgp
       Tn=Tn+dTn

       Tabc=(((ep*ep)/(1.0d-18))**(0.25))*0.246*Tg
       Tabc=Tabc*(((1.0/Tg)-(1.0/Tinitial))**(0.25))

       if(Tg.lt.10) then

c start calculation of lam1

       lam1=0
       da1=Tn*14.0/1000.0
       a1=0
       do 7 i7=1,1000,+1
       Ee=Q+a1
       dlam1=da1*Ee*Ee*a1*a1*Awk
       dlam1=dlam1/(1.0+dexp(a1/Tn))
       dlam1=dlam1/(1.0+dexp(-Ee/Tg))
       lam1=lam1+dlam1
       a1=a1+da1
7      continue

c end calculation of lam1
c start calculation of lam2

       lam2=0
       da2=Tg*14.0/1000.0
       a2=0
       do 8 i8=1,1000,+1
       Ee=((a2*a2+me*me)**(0.5))
       En=Q+Ee
       dlam2=da2*En*En*a2*a2*Awk
       dlam2=dlam2/(1.0+dexp(Ee/Tg))
       dlam2=dlam2/(1.0+dexp(-En/Tn))
       lam2=lam2+dlam2
       a2=a2+da2
8      continue

c end calculation of lam2
c start calculation of lam3

       lam3=0
       da3=Tg*14.0/1000.0
       a3=ext3+da3
       do 9 i9=1,1000,+1
       Ee=((a3*a3+me*me)**(0.5))
       En=Ee-Q
       dlam3=da3*En*En*a3*a3*Awk
       dlam3=dlam3/(1.0+dexp(Ee/Tg))
       dlam3=dlam3/(1.0+dexp(-En/Tn))
       lam3=lam3+dlam3
       a3=a3+da3
9      continue

c end calculation of lam3
c start calculation of lam4

       lam4=0
       da4=Tn*14.0/1000.0
       a4=ext4+da4
       do 10 i10=1,1000,+1
       Ee=a4-Q
       dlam4=da4*Ee*Ee*a4*a4*Awk
       dlam4=dlam4/(1.0+dexp(a4/Tn))
       dlam4=dlam4/(1.0+dexp(-Ee/Tg))
       lam4=lam4+dlam4
       a4=a4+da4
10     continue

       dXn=((-(lam1+lam2)*Xn)+((lam3+lam4)*(1.0-Xn)))*dt
       Xn=Xn+dXn

       else 
       goto 92
       endif

c end calculation of lam4

92     if(i.gt.100) then
c       write(6,*) 't=', t*(6.582d-22)
c       write(6,*) 'Tg=', Tg, 'Tn=', Tn, 'Tgp=', Tgp, 'Tgpp=', Tabc
c       write(6,*) 'Xn=', Xn
c       write(6,*) 'ratio=', Tgp/Tg
       i=1
       endif
       i=i+1
c set time step according to neutrino temperature

       if(Tn.gt.0.6) then
       dt=40*(((10.0/Tn)*(10.0/Tn)*(1.5E21))/(1280000.0))
       else
       dt=25*(((10.0/Tn)*(10.0/Tn)*(1.5E21))/(2510000.0))
       endif
       t=t+dt
       if(Tg.lt.0.07)then
       goto 11
       endif
6      continue
11     tp=t*6.582d-22 
c       write(6,*) 'Yp=', 2*Xn*exp(-tp/886.8), 't=', tp
       Yp=2*Xn*exp(-tp/886.8)
c       if(j.eq.1) then
c       Ype0=Yp
c       endif
       write(6,*) ep, (Yp-0.24774297)/0.013
       ep=ep+depsilon
93     continue
       end
