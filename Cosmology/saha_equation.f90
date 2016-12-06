       Implicit none
       real*8 y,A,yp,me,mf1,mf2,x,ion,ep,ion1,alphap,T,Tp
       integer i

       me=511000.0 !eV
       mf1=10000000.0 !eV
       mf2=2.0     !GeV
       alphap=1/(10.0)
       ion=11.0 !eV
       ion1=ion*((alphap*137.0)**(2.0))*(mf1/me)
       ep=4.4d-2
       if(mf1.gt.me) then
       x=mf1
       else
       x=me
       endif
       A=3.5d-7
       A=A*(((1.0d-9)/ep)**(1.5))
       A=A*((x/me)**(0.75))
       A=A/mf2
       A=A*((ion1/mf1)**(1.5))
       yp=1
       do 1 i=1,3,+1
       y=((1.5)*(dlog(yp)))+dlog(90.0/A)
       yp=y
1      continue

       Tp=ion1/y
       T=Tp/(0.31*((ep/(1.0d-9))**(0.5))*((me/x)**(0.25)))
       write(6,*) 'y=',y,'Tp=',Tp
c       write(6,*) 'T=',ion/(y*0.31*((ep/1.0d-9)**(0.5))*((me/x)**(0.25)))
       write(6,*) 'T=',T


       end

