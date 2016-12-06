       Implicit none
       real*8 nr(0:6371),Tr(0:6371),rhor(0:6371),gr(0:6371)
       real*8 dr,pi,mf2,mp,n0,gnp,dTr,dnr,dgr
       real*8 mf2k,conv,drp,temp,jgr,vrot,vrot1,conv1
       real*8 Zp,alphap,mf1,conv2,mf2cm,ext1,ext2,re,q,dq
       real*8 lhs,dlhs,rpp,r0,r0p,insl,rhs,conv3
       integer i1,i2,i3,i4,i5,i6,irpp

       mf2=40.0 !GeV
       pi=3.1415927
       n0=1.0d14 !cm^-3
       dr=1.0000 !km
       gnp=6.67d-5 !Netwon's constant adjusted for different units
       conv=1.29d-4 !converts mf2(GeV) to mf2(g)/K_B
       conv1=3.0d5 !speed of light in km/s
       conv2=5.06d13 !converts mf2(GeV) to mf2(cm^-1)
       conv3=2.56d27 !converts GeV^2 to cm^-2
       drp=dr*1000.0 !m
       vrot=220 !km/s
       Zp=60.0
       alphap=1.0d-2
       mf1=1.0d-3 !GeV
       re=6371.0 !km
c       r0=4500.0 !km

c convert quantities and set the initial conditions
       
       mf2k=mf2*conv
       mf2cm=mf2*conv2
       vrot1=vrot/conv1
       gr(0)=0
       nr(0)=n0
       
c temperature and density profile

       do 1 i1=0,6371,+1
       if(i1.ge.6347) then
       Tr(i1)=398470.9354-62.4974*i1
       rhor(i1)=965.939165-0.151615*i1
       else if (i1.ge.3630) then
       Tr(i1)=5730-0.6192*i1
       rhor(i1)=8.2086-0.00072*i1
       else if (i1.ge.3480) then
       Tr(i1)=5730-0.6192*i1
       rhor(i1)=117.036-0.0307*i1
       else if (i1.ge.1221) then
       Tr(i1)=5730-0.6192*i1
       rhor(i1)=14.028-0.0011*i1
       else
       Tr(i1)=5730-0.6192*i1
       rhor(i1)=13.0885
       endif
1      continue

c start calculation of g(r), gravitational acceleration inside earth

       jgr=0.0
       do 2 i2=1,6371,+1
       jgr=jgr+4*pi*i2*i2*rhor(i2)*dr
       gr(i2)=jgr*gnp/(i2*i2)
2      continue      

c start calculation of n(r)

       do 3 i3=1,6370,+1
       if(i3.ge.6347) then
       dTr=62.4974
       else
       dTr=0.6192
       endif
       temp=(mf2k*gr(i3-1)*drp)+dTr
       dnr=-nr(i3-1)/Tr(i3-1)
       dnr=dnr*temp
       nr(i3)=nr(i3-1)+dnr
c       write(6,*) i3, nr(i3)
3      continue

c loop on distance of closest approach

       do 5 i5=0,6371,+1
       r0=i5*1.0
       ext1=-((re*re-r0*r0)**(0.5))
       ext2=((re*re-r0*r0)**(0.5))
       r0p=r0*1.0d5 !cm
       dq=((ext2-ext1)*1.0d5)/6371.0 !cm
       q=ext1*1.0d5 !cm       
       lhs=0.0

c       do 6 i6=1,2000,+1
c       rpp=(r0*r0+q*q)**(0.5)
c       irpp=nint(rpp)
c       q=q+dq
c       write(6,*) irpp
c6      continue

c calculate integral of n(x) in dx

       do 4 i4=0,6371,+1
       rpp=(r0p*r0p+q*q)**(0.5)
       rpp=rpp/(1.0d5)
       irpp=nint(rpp)
       dlhs=nr(irpp)*dq
       lhs=lhs+dlhs
       q=q+dq
4      continue

c calculate rhs to compare with integral of n(x) dx

       insl=(mf2*vrot1)/(mf1*alphap)
       rhs=mf2*mf2*vrot1*vrot1*vrot1*vrot1
       rhs=rhs/(8.0*pi*Zp*Zp*Zp*Zp)
       rhs=rhs/(alphap*alphap)
       rhs=rhs/(log(insl))
       rhs=rhs*conv3
       
       write(6,*) i5, lhs/(1.0d13)-rhs/(1.0d13)
5      continue
       end
       
       
