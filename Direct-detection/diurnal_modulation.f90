       Implicit none
       real*8 mf2,alphap,Zp,c1,c2,c3,r0,re,dmin2,Er,mt
       real*8 f,g,halo,vrot,ve,lat,pi
       real*8 psi,phi,dphi,cth,dcth,v,dv,vmin,cpsi,theta
       real*8 ipsi,modul,v0,v02,mbar,fa,t,dipsi,dmodul,dt
       integer north,i1,i2,i3,i4

       t=0.0 !time in sidereal days units
       dt=0.01

       do 1 i1=1,101,+1 !time loop
       mf2=10.0 !GeV
       alphap=1.0d-2
       Zp=10.0
       Er=2.0d-6 !GeV
       mt=22.98*0.94 !GeV
       pi=3.1415927
       lat=37.03*pi/180.0 !latitude of detector in radians
       vrot=220.0 !km/s
       ve=vrot+12.0 !km/s
       halo=43.0*pi/180.0 !angle earth motion wrt spin rad
       north=0 !1=northern hemisphere,0=southern hemisphere
       re=6371.0 !km
       
       mbar=mf2/Zp !mean mass
       v0=vrot*((mbar/mf2)**(0.5)) !velocity dispersion
       v02=v0*v0
       c1=mf2/10.0
       c2=alphap/(1.0d-3)
       c3=Zp/10.0
       r0=5000.0*(c1**(-0.55))*(c2**(0.06))*(c3**(0.15)) !km
       vmin=(((mt+mf2)*(mt+mf2)*Er*0.5)/(mt*mf2*mf2))**(0.5)
       vmin=vmin*3.0d5 !km/s
       
       dphi=2.0*pi/200.0
       dcth=2.0/200.0
       dv=500.0/200.0
       v=vmin

       if(north.eq.1) then
       cpsi=dcos(lat)*dsin(2*pi*(t))*dsin(halo)+dsin(lat)*dcos(halo)
       else
       cpsi=dcos(lat)*dsin(2*pi*(t))*dsin(halo)-dsin(lat)*dcos(halo)
       endif

       psi=acos(cpsi)

       ipsi=0.0
       modul=0.0
       
       do 2 i2=1,200,+1 !v loop
       cth=-1.0
       do 3 i3=1,200,+1 !cos theta loop
       phi=0.0
       do 4 i4=1,200,+1 !phi loop
       theta=acos(cth)
       f=dsin(theta)*dsin(psi)*dsin(phi)-cth*dcos(psi)

       if(f.gt.0) then
       dmin2=(1-(f*f))*re*re
       else
       dmin2=re*re
       endif
       
       if(dmin2.lt.(r0*r0)) then
       g=0.0
       else
       g=1.0
       endif

       dipsi=(dexp(-v*v/v02))
       dipsi=dipsi*dexp(-ve*ve/v02)
       dipsi=dipsi*dexp(2.0*v*ve*cth/v02)
       dipsi=dipsi*v*dphi*dcth*dv
       dmodul=dipsi*g
       ipsi=ipsi+dipsi
       modul=modul+dmodul
       phi=phi+dphi
4      continue
       cth=cth+dcth
3      continue
       v=v+dv
2      continue
       write(6,*) t,100*(1-modul/ipsi)
       t=t+dt
1      continue
       end
