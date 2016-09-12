CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCC Numerator functions CCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm2num1e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin,s34spinc
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag,s34spinc
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      taux1 = 512.D0/3.D0*s3l1*s4l2+512.D0/3.D0*s3l2*s4l1
      taux2 = 2048.D0/3.D0*s3l1*s4l2-2048.D0/3.D0*s3l2*s4l1

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = taux2

CCC Jet functions

      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet1 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t1 = s14f(z,l1,l2,l3,l4)
      t2 = s24f(z,l1,l2,l3,l4)
      t3 = s13f(z,l1,l2,l3,l4)
      t4 = s23f(z,l1,l2,l3,l4)
      tFjet2 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet3 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = s14f(z,l1,l2,l3,l4)
      t2 = s24f(z,l1,l2,l3,l4)
      t3 = s13f(z,l1,l2,l3,l4)
      t4 = s23f(z,l1,l2,l3,l4)
      tFjet4 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff4)

CCC Sum

      rm2num1e0 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff
     #4
      RETURN
      END

      double precision function rm2num1e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin,s34spinc
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag,s34spinc
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = s13f(z,l1,l2,l3,l4)
      t4 = s14f(z,l1,l2,l3,l4)
      t7 = s23f(z,l1,l2,l3,l4)
      t10 = s24f(z,l1,l2,l3,l4)
      taux1 = 512.D0/3.D0*z*t1+512.D0/3.D0*z*t4+512.D0/3.D0*z*t7+512.D0/
     #3.D0*z*t10-1024.D0/3.D0*s3l1*s4l2-1024.D0/3.D0*s3l2*s4l1-512.D0/3.
     #D0*z*(-1+z)
      taux2 = -10240.D0/3.D0*s3l1*s4l2+10240.D0/3.D0*s3l2*s4l1

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = taux2

CCC Jet functions

      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet1 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t1 = s14f(z,l1,l2,l3,l4)
      t2 = s24f(z,l1,l2,l3,l4)
      t3 = s13f(z,l1,l2,l3,l4)
      t4 = s23f(z,l1,l2,l3,l4)
      tFjet2 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet3 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = s14f(z,l1,l2,l3,l4)
      t2 = s24f(z,l1,l2,l3,l4)
      t3 = s13f(z,l1,l2,l3,l4)
      t4 = s23f(z,l1,l2,l3,l4)
      tFjet4 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff4)

CCC Sum

      rm2num1e1 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff
     #4
      RETURN
      END

      double precision function rm2num1e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin,s34spinc
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag,s34spinc
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = s13f(z,l1,l2,l3,l4)
      t4 = s14f(z,l1,l2,l3,l4)
      t7 = s23f(z,l1,l2,l3,l4)
      t10 = s24f(z,l1,l2,l3,l4)
      taux1 = -1024.D0/3.D0*z*t1-1024.D0/3.D0*z*t4-1024.D0/3.D0*z*t7-102
     #4.D0/3.D0*z*t10+512.D0/3.D0*s3l1*s4l2+512.D0/3.D0*s3l2*s4l1+1024.D
     #0/3.D0*z*(-1+z)
      taux2 = 6144*s3l1*s4l2-6144*s3l2*s4l1

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = taux2

CCC Jet functions

      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet1 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t1 = s14f(z,l1,l2,l3,l4)
      t2 = s24f(z,l1,l2,l3,l4)
      t3 = s13f(z,l1,l2,l3,l4)
      t4 = s23f(z,l1,l2,l3,l4)
      tFjet2 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet3 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = s14f(z,l1,l2,l3,l4)
      t2 = s24f(z,l1,l2,l3,l4)
      t3 = s13f(z,l1,l2,l3,l4)
      t4 = s23f(z,l1,l2,l3,l4)
      tFjet4 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff4)

CCC Sum

      rm2num1e2 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff
     #4
      RETURN
      END

      double precision function rm2num1e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin,s34spinc
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag,s34spinc
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = s13f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t5 = s23f(z,l1,l2,l3,l4)
      t7 = s24f(z,l1,l2,l3,l4)
      taux1 = 512.D0/3.D0*z*t1+512.D0/3.D0*z*t3+512.D0/3.D0*z*t5+512.D0/
     #3.D0*z*t7-512.D0/3.D0*z*(-1+z)
      taux2 = -14336.D0/3.D0*s3l1*s4l2+14336.D0/3.D0*s3l2*s4l1

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = taux2

CCC Jet functions

      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet1 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t1 = s14f(z,l1,l2,l3,l4)
      t2 = s24f(z,l1,l2,l3,l4)
      t3 = s13f(z,l1,l2,l3,l4)
      t4 = s23f(z,l1,l2,l3,l4)
      tFjet2 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet3 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = s14f(z,l1,l2,l3,l4)
      t2 = s24f(z,l1,l2,l3,l4)
      t3 = s13f(z,l1,l2,l3,l4)
      t4 = s23f(z,l1,l2,l3,l4)
      tFjet4 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff4)

CCC Sum

      rm2num1e3 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff
     #4
      RETURN
      END

      double precision function rm2num1e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin,s34spinc
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag,s34spinc
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      taux1 = 0
      taux2 = 4096.D0/3.D0*s3l1*s4l2-4096.D0/3.D0*s3l2*s4l1

CCC Coefficients

      tcoff1 = taux1
      tcoff2 = taux1
      tcoff3 = taux2
      tcoff4 = taux2

CCC Jet functions

      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet1 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff1)
      t1 = s14f(z,l1,l2,l3,l4)
      t2 = s24f(z,l1,l2,l3,l4)
      t3 = s13f(z,l1,l2,l3,l4)
      t4 = s23f(z,l1,l2,l3,l4)
      tFjet2 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff2)
      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet3 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff3)
      t1 = s14f(z,l1,l2,l3,l4)
      t2 = s24f(z,l1,l2,l3,l4)
      t3 = s13f(z,l1,l2,l3,l4)
      t4 = s23f(z,l1,l2,l3,l4)
      tFjet4 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s4l1,s3l1,s1l2,s2l2,s4l2,s3
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff4)

CCC Sum

      rm2num1e4 = tFjet1*tcoff1+tFjet2*tcoff2+tFjet3*tcoff3+tFjet4*tcoff
     #4
      RETURN
      END

      double precision function rm2num2e0(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin,s34spinc
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag,s34spinc
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = 128.D0/3.D0*s1l2
      t2 = 64.D0/3.D0*s2l2
      t3 = 64.D0/3.D0*s4l2
      t9 = s24f(z,l1,l2,l3,l4)
      t15 = 64.D0/3.D0*s2l2*s4l1
      t16 = 128.D0/3.D0*z
      t18 = s13f(z,l1,l2,l3,l4)
      t20 = 64.D0/3.D0*s3l2
      t26 = s23f(z,l1,l2,l3,l4)
      t32 = 64.D0/3.D0*s2l2*s3l1
      t34 = s14f(z,l1,l2,l3,l4)
      t36 = -1+z
      t37 = 128.D0/3.D0*t36
      t61 = 128.D0/3.D0*s2l2
      taux3 = ((-t1+t2-t3)*s1l1+(64.D0/3.D0*s2l1-64.D0/3.D0*s4l1)*s1l2-1
     #28.D0/3.D0*z*t9+64.D0/3.D0*s2l1*s4l2+t15-t16)*t18+((-t1+t2-t20)*s1
     #l1+(64.D0/3.D0*s2l1-64.D0/3.D0*s3l1)*s1l2-128.D0/3.D0*z*t26+64.D0/
     #3.D0*s2l1*s3l2+t32-t16)*t34+(t37*s1l2+(64.D0/3.D0*s2l2+64.D0/3.D0*
     #s4l2)*t26+(64.D0/3.D0*s2l2+64.D0/3.D0*s3l2)*t9-128.D0/3.D0*z*s2l2-
     #t20-t3)*s1l1+((64.D0/3.D0*s2l1+64.D0/3.D0*s4l1)*t26+(64.D0/3.D0*s2
     #l1+64.D0/3.D0*s3l1)*t9-128.D0/3.D0*z*s2l1-64.D0/3.D0*s3l1-64.D0/3.
     #D0*s4l1)*s1l2+((-t61-t3)*s2l1-t15-t16)*t26+((-t61-t20)*s2l1-t32-t1
     #6)*t9+(t37*s2l2-t20-t3)*s2l1+(-64.D0/3.D0*s3l1-64.D0/3.D0*s4l1)*s2
     #l2-128.D0/3.D0*s3l1*s4l2-128.D0/3.D0*s3l2*s4l1+128.D0/3.D0*z*t36
      t2 = -256.D0/3.D0*s2l2-256.D0/3.D0*s4l2
      t5 = 256.D0/3.D0*s2l1+256.D0/3.D0*s4l1
      t7 = s2l1*s4l2
      t9 = s2l2*s4l1
      t12 = s13f(z,l1,l2,l3,l4)
      t15 = -256.D0/3.D0*s2l2-256.D0/3.D0*s3l2
      t18 = 256.D0/3.D0*s2l1+256.D0/3.D0*s3l1
      t20 = s2l1*s3l2
      t22 = s2l2*s3l1
      t25 = s14f(z,l1,l2,l3,l4)
      t27 = s23f(z,l1,l2,l3,l4)
      t29 = s24f(z,l1,l2,l3,l4)
      taux4 = (t2*s1l1+t5*s1l2+256.D0/3.D0*t7-256.D0/3.D0*t9)*t12+(t15*s
     #1l1+t18*s1l2+256.D0/3.D0*t20-256.D0/3.D0*t22)*t25+(t2*t27+t15*t29-
     #512.D0/3.D0*s2l2-256.D0/3.D0*s3l2-256.D0/3.D0*s4l2)*s1l1+(t5*t27+t
     #18*t29+512.D0/3.D0*s2l1+256.D0/3.D0*s3l1+256.D0/3.D0*s4l1)*s1l2+(2
     #56.D0/3.D0*t7-256.D0/3.D0*t9)*t27+(256.D0/3.D0*t20-256.D0/3.D0*t22
     #)*t29+(256.D0/3.D0*s3l2+256.D0/3.D0*s4l2)*s2l1+(-256.D0/3.D0*s3l1-
     #256.D0/3.D0*s4l1)*s2l2

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux4

CCC Jet functions

      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet5 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet6 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff6)

CCC Sum

      rm2num2e0 = tFjet5*tcoff5+tFjet6*tcoff6
      RETURN
      END

      double precision function rm2num2e1(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin,s34spinc
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag,s34spinc
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = s1l1*s2l2
      t2 = s1l2*s2l1
      t3 = s24f(z,l1,l2,l3,l4)
      t4 = z*t3
      t5 = s2l1*s4l2
      t6 = s2l2*s4l1
      t9 = s13f(z,l1,l2,l3,l4)
      t10 = t9**2
      t14 = s23f(z,l1,l2,l3,l4)
      t15 = z*t14
      t17 = 64.D0/3.D0*t4
      t18 = s3l2+s4l2
      t21 = s3l1+s4l1
      t25 = s14f(z,l1,l2,l3,l4)
      t27 = 128.D0/3.D0*s2l2
      t28 = 64.D0/3.D0*s4l2
      t33 = 1-z
      t34 = 128.D0/3.D0*t33
      t35 = t34*s2l2
      t36 = 64.D0/3.D0*s3l2
      t39 = 128.D0/3.D0*s2l1
      t40 = 64.D0/3.D0*s4l1
      t45 = t34*s2l1
      t46 = 64.D0/3.D0*s3l1
      t52 = 128.D0/3.D0*s4l2*s4l1
      t53 = 128.D0/3.D0*z
      t56 = t3**2
      t58 = 64.D0/3.D0*z*t56
      t59 = s3l1*s4l2
      t60 = 128.D0/3.D0*t59
      t61 = s3l2*s4l1
      t62 = 128.D0/3.D0*t61
      t65 = 64.D0/3.D0*z*(-3+z)
      t68 = z*s4l2
      t71 = z*s4l1
      t74 = 64*t59
      t75 = 64*t61
      t76 = -t33
      t79 = 128.D0/3.D0*z*t76
      t80 = (128.D0/3.D0*t1+128.D0/3.D0*t2-64.D0/3.D0*t15-t17+128.D0/3.D
     #0*t18*s2l1+128.D0/3.D0*t21*s2l2)*t25+((t27+t28)*t14+128.D0/3.D0*s2
     #l2*t3+t35-t36-t28)*s1l1+((t39+t40)*t14+128.D0/3.D0*s2l1*t3+t45-t46
     #-t40)*s1l2+(-t17+64.D0/3.D0*t5+64.D0/3.D0*t6+t52-t53)*t14-t58+(-t6
     #0-t62+t65)*t3-64.D0/3.D0*t68*s2l1-64.D0/3.D0*t71*s2l2-t74-t75-t52+
     #t79
      t82 = s2l1*s3l2
      t83 = s2l2*s3l1
      t86 = t25**2
      t100 = t14**2
      t102 = 64.D0/3.D0*z*t100
      t111 = z*s3l2
      t114 = z*s3l1
      t117 = 128.D0/3.D0*s3l2
      t120 = (-t117-64*s4l2)*s3l1
      t138 = z**2
      t140 = 64.D0/3.D0-t53+64.D0/3.D0*t138
      t142 = 64.D0/3.D0*t76
      t143 = t142*s3l2
      t144 = t142*s4l2
      t163 = t142*s3l1
      t164 = t142*s4l1
      t169 = -64.D0/3.D0*t18*s2l1
      t172 = -64.D0/3.D0*t21*s2l2
      t181 = 64*t76
      s2 = (64.D0/3.D0*t1+64.D0/3.D0*t2-64.D0/3.D0*t4+64.D0/3.D0*t5+64.D
     #0/3.D0*t6-64.D0/3.D0*z)*t10+t80*t9+(64.D0/3.D0*t1+64.D0/3.D0*t2-64
     #.D0/3.D0*t15+64.D0/3.D0*t82+64.D0/3.D0*t83-64.D0/3.D0*z)*t86+((128
     #.D0/3.D0*s2l2*t14+(t27+t36)*t3+t35-t36-t28)*s1l1+(128.D0/3.D0*s2l1
     #*t14+(t39+t46)*t3+t45-t46-t40)*s1l2-t102+(-t17-t60-t62+t65)*t14+(6
     #4.D0/3.D0*t82+64.D0/3.D0*t83+128.D0/3.D0*s3l2*s3l1-t53)*t3-64.D0/3
     #.D0*t111*s2l1-64.D0/3.D0*t114*s2l2+t120-t75+t79)*t25
      s1 = s2+((64.D0/3.D0*s2l2+64.D0/3.D0*s4l2)*t100+((128.D0/3.D0*s2l2
     #+128.D0/3.D0*s3l2+128.D0/3.D0*s4l2)*t3+t35-64.D0/3.D0*t68)*t14+(64
     #.D0/3.D0*s2l2+64.D0/3.D0*s3l2)*t56+(t35-64.D0/3.D0*t111)*t3+t140*s
     #2l2+t143+t144)*s1l1+((64.D0/3.D0*s2l1+64.D0/3.D0*s4l1)*t100+((128.
     #D0/3.D0*s2l1+128.D0/3.D0*s3l1+128.D0/3.D0*s4l1)*t3+t45-64.D0/3.D0*
     #t71)*t14+(64.D0/3.D0*s2l1+64.D0/3.D0*s3l1)*t56+(t45-64.D0/3.D0*t11
     #4)*t3+t140*s2l1+t163+t164)*s1l2-t102+(t169+t172-t74-t75-t52+t79)*t
     #14
      taux3 = s1-t58+(t169+t172+t120-t75+t79)*t3+(t143+t144)*s2l1+(t163+
     #t164)*s2l2+(-t117+t181*s4l2)*s3l1+t181*s4l1*s3l2-t52-64.D0/3.D0*z*
     #(1-4*z+t138)
      t1 = s1l1*s2l2
      t2 = s1l2*s2l1
      t3 = s2l1*s4l2
      t4 = s2l2*s4l1
      t7 = s13f(z,l1,l2,l3,l4)
      t8 = t7**2
      t19 = s14f(z,l1,l2,l3,l4)
      t21 = 512.D0/3.D0*s2l2
      t22 = 256.D0/3.D0*s4l2
      t24 = s23f(z,l1,l2,l3,l4)
      t26 = 512.D0/3.D0*s3l2
      t28 = s24f(z,l1,l2,l3,l4)
      t31 = 256*z*s2l2
      t32 = 256.D0/3.D0*s3l2
      t33 = 512.D0/3.D0*s4l2
      t36 = 512.D0/3.D0*s2l1
      t37 = 256.D0/3.D0*s4l1
      t40 = 512.D0/3.D0*s3l1
      t44 = 256*z*s2l1
      t45 = 256.D0/3.D0*s3l1
      t46 = 512.D0/3.D0*s4l1
      t56 = s3l1*s4l2
      t57 = 512.D0/3.D0*t56
      t58 = s3l2*s4l1
      t59 = 512.D0/3.D0*t58
      t63 = 512.D0/3.D0*z*s4l2
      t67 = 512.D0/3.D0*z*s4l1
      t70 = 256.D0/3.D0*t56
      t71 = 256.D0/3.D0*t58
      t74 = s2l1*s3l2
      t75 = s2l2*s3l1
      t78 = t19**2
      t102 = 512.D0/3.D0*z*s3l2
      t106 = 512.D0/3.D0*z*s3l1
      t113 = t24**2
      t123 = t28**2
      t128 = z**2
      t130 = 4*z-2*t128+4
      t133 = -z-2
      t134 = t133*s3l2
      t136 = t133*s4l2
      t156 = t133*s3l1
      t158 = t133*s4l1
      s2 = (-256.D0/3.D0*t1+256.D0/3.D0*t2+256.D0/3.D0*t3-256.D0/3.D0*t4
     #)*t8+((-512.D0/3.D0*t1+512.D0/3.D0*t2+(256*s3l2+256*s4l2)*s2l1+(-2
     #56*s3l1-256*s4l1)*s2l2)*t19+((-t21-t22)*t24+(-t21-t26-t22)*t28+t31
     #-t32+t33)*s1l1+((t36+t37)*t24+(t36+t40+t37)*t28-t44+t45-t46)*s1l2+
     #(256.D0/3.D0*t3-256.D0/3.D0*t4)*t24+((t32+t33)*s2l1+(-t45-t46)*s2l
     #2+t57-t59)*t28+(t32-t63)*s2l1+(-t45+t67)*s2l2+t70-t71)*t7
      s1 = s2+(-256.D0/3.D0*t1+256.D0/3.D0*t2+256.D0/3.D0*t74-256.D0/3.D
     #0*t75)*t78+(((-t21-t32-t33)*t24+(-t21-t32)*t28+t31+t26-t22)*s1l1+(
     #(t36+t45+t46)*t24+(t36+t45)*t28-t44-t40+t37)*s1l2+((t26+t22)*s2l1+
     #(-t40-t37)*s2l2-t57+t59)*t24+(256.D0/3.D0*t74-256.D0/3.D0*t75)*t28
     #+(-t102+t22)*s2l1+(t106-t37)*s2l2-t70+t71)*t19+((-256.D0/3.D0*s2l2
     #-256.D0/3.D0*s4l2)*t113+((-t21-256*s3l2-256*s4l2)*t28+t31-t32+t63)
     #*t24+(-256.D0/3.D0*s2l2-256.D0/3.D0*s3l2)*t123+(t31+t102-t22)*t28+
     #256.D0/3.D0*t130*s2l2-256.D0/3.D0*t134-256.D0/3.D0*t136)*s1l1
      taux4 = s1+((256.D0/3.D0*s2l1+256.D0/3.D0*s4l1)*t113+((t36+256*s3l
     #1+256*s4l1)*t28-t44+t45-t67)*t24+(256.D0/3.D0*s2l1+256.D0/3.D0*s3l
     #1)*t123+(-t44-t106+t37)*t28-256.D0/3.D0*t130*s2l1+256.D0/3.D0*t156
     #+256.D0/3.D0*t158)*s1l2+((t32-t33)*s2l1+(-t45+t46)*s2l2-t70+t71)*t
     #24+((-t26+t22)*s2l1+(t40-t37)*s2l2+t70-t71)*t28+(256.D0/3.D0*t134+
     #256.D0/3.D0*t136)*s2l1+(-256.D0/3.D0*t156-256.D0/3.D0*t158)*s2l2

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux4

CCC Jet functions

      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet5 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet6 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff6)

CCC Sum

      rm2num2e1 = tFjet5*tcoff5+tFjet6*tcoff6
      RETURN
      END

      double precision function rm2num2e2(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin,s34spinc
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag,s34spinc
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = s1l1*s2l2
      t2 = 64.D0/3.D0*t1
      t3 = s1l2*s2l1
      t4 = 64.D0/3.D0*t3
      t5 = s24f(z,l1,l2,l3,l4)
      t6 = z*t5
      t8 = s2l1*s4l2
      t9 = 64.D0/3.D0*t8
      t10 = s2l2*s4l1
      t11 = 64.D0/3.D0*t10
      t12 = 64*z
      t14 = s13f(z,l1,l2,l3,l4)
      t15 = t14**2
      t19 = s23f(z,l1,l2,l3,l4)
      t20 = z*t19
      t22 = 320.D0/3.D0*t6
      t31 = s14f(z,l1,l2,l3,l4)
      t33 = 128.D0/3.D0*s2l2
      t34 = 64.D0/3.D0*s4l2
      t37 = 64.D0/3.D0*s3l2
      t40 = 1-z
      t41 = 128.D0/3.D0*t40
      t42 = t41*s2l2
      t44 = ((t33+t34)*t19+(t33+t37)*t5+t42)*s1l1
      t45 = 128.D0/3.D0*s2l1
      t46 = 64.D0/3.D0*s4l1
      t49 = 64.D0/3.D0*s3l1
      t52 = t41*s2l1
      t54 = ((t45+t46)*t19+(t45+t49)*t5+t52)*s1l2
      t57 = 512.D0/3.D0*z
      t60 = t5**2
      t62 = 64*z*t60
      t65 = z*(-10+4*z)
      t69 = 64.D0/3.D0*t40
      t70 = t69*s4l2
      t72 = t69*s4l1
      t74 = s3l1*s4l2
      t75 = 64.D0/3.D0*t74
      t76 = s3l2*s4l1
      t77 = 64.D0/3.D0*t76
      t80 = -7*z*t40
      t81 = 64.D0/3.D0*t80
      t82 = (128.D0/3.D0*t1+128.D0/3.D0*t3-320.D0/3.D0*t20-t22+(64.D0/3.
     #D0*s3l2+64.D0/3.D0*s4l2)*s2l1+(64.D0/3.D0*s3l1+64.D0/3.D0*s4l1)*s2
     #l2-128*z)*t31+t44+t54+(-t22+t9+t11+128.D0/3.D0*s4l2*s4l1-t57)*t19-
     #t62+(64.D0/3.D0*t8+64.D0/3.D0*t10+64.D0/3.D0*t65)*t5+t70*s2l1+t72*
     #s2l2-t75-t77+t81
      t86 = 64.D0/3.D0*s2l1*s3l2
      t88 = 64.D0/3.D0*s2l2*s3l1
      t90 = t31**2
      t92 = t19**2
      t94 = 64*z*t92
      t102 = t69*s3l2
      t104 = t69*s3l1
      t121 = z**2
      t123 = 64.D0/3.D0-128.D0/3.D0*z+64.D0/3.D0*t121
      t148 = -t69
      taux3 = (t2+t4-64*t6+t9+t11-t12)*t15+t82*t14+(t2+t4-64*t20+t86+t88
     #-t12)*t90+(t44+t54-t94+(-t22+t86+t88+64.D0/3.D0*t65)*t19+(t86+t88+
     #128.D0/3.D0*s3l2*s3l1-t57)*t5+t102*s2l1+t104*s2l2-t75-t77+t81)*t31
     #+((64.D0/3.D0*s2l2+64.D0/3.D0*s4l2)*t92+((t33+t37+t34)*t5+t42+t70)
     #*t19+(64.D0/3.D0*s2l2+64.D0/3.D0*s3l2)*t60+(t42+t102)*t5+t123*s2l2
     #)*s1l1+((64.D0/3.D0*s2l1+64.D0/3.D0*s4l1)*t92+((t45+t49+t46)*t5+t5
     #2+t72)*t19+(64.D0/3.D0*s2l1+64.D0/3.D0*s3l1)*t60+(t52+t104)*t5+t12
     #3*s2l1)*s1l2-t94+(-128*t6-t75-t77+t81)*t19-t62+(-64.D0/3.D0*t74-64
     #.D0/3.D0*t76+64.D0/3.D0*t80)*t5+t148*s4l2*s3l1+t148*s4l1*s3l2-64.D
     #0/3.D0*z*(4-8*z+4*t121)
      t1 = s1l1*s2l2
      t2 = s1l2*s2l1
      t3 = s2l1*s4l2
      t4 = s2l2*s4l1
      t7 = s13f(z,l1,l2,l3,l4)
      t8 = t7**2
      t12 = -s3l2-s4l2
      t15 = s3l1+s4l1
      t19 = s14f(z,l1,l2,l3,l4)
      t21 = 512.D0/3.D0*s2l2
      t22 = 256.D0/3.D0*s4l2
      t24 = s23f(z,l1,l2,l3,l4)
      t26 = 256.D0/3.D0*s3l2
      t27 = -t21+t26+t22
      t28 = s24f(z,l1,l2,l3,l4)
      t33 = 512.D0/3.D0*s2l1
      t34 = 256.D0/3.D0*s4l1
      t37 = 256.D0/3.D0*s3l1
      t38 = t33-t37-t34
      t47 = 256.D0/3.D0*t12*s2l1
      t49 = 256.D0/3.D0*t15*s2l2
      t50 = s3l1*s4l2
      t51 = 256.D0/3.D0*t50
      t52 = s3l2*s4l1
      t53 = 256.D0/3.D0*t52
      t56 = -1+z
      t57 = t56*s4l2
      t61 = t56*s4l1
      t65 = 512.D0/3.D0*t50
      t66 = 512.D0/3.D0*t52
      t69 = s2l1*s3l2
      t70 = s2l2*s3l1
      t73 = t19**2
      t92 = t56*s3l2
      t96 = t56*s3l1
      t104 = t24**2
      t114 = t28**2
      t120 = z**2
      t121 = -2*z-3+t120
      t124 = z*s3l2
      t126 = z*s4l2
      t147 = z*s3l1
      t149 = z*s4l1
      s2 = (-256.D0/3.D0*t1+256.D0/3.D0*t2+256.D0/3.D0*t3-256.D0/3.D0*t4
     #)*t8+((-512.D0/3.D0*t1+512.D0/3.D0*t2+512.D0/3.D0*t12*s2l1+512.D0/
     #3.D0*t15*s2l2)*t19+((-t21-t22)*t24+t27*t28-t21+512.D0/3.D0*s3l2)*s
     #1l1+((t33+t34)*t24+t38*t28+t33-512.D0/3.D0*s3l1)*s1l2+(256.D0/3.D0
     #*t3-256.D0/3.D0*t4)*t24+(t47+t49-t51+t53)*t28+(-256.D0/3.D0*s3l2+2
     #56.D0/3.D0*t57)*s2l1+(256.D0/3.D0*s3l1-256.D0/3.D0*t61)*s2l2-t65+t
     #66)*t7
      s1 = s2+(-256.D0/3.D0*t1+256.D0/3.D0*t2+256.D0/3.D0*t69-256.D0/3.D
     #0*t70)*t73+((t27*t24+(-t21-t26)*t28-t21+512.D0/3.D0*s4l2)*s1l1+(t3
     #8*t24+(t33+t37)*t28+t33-512.D0/3.D0*s4l1)*s1l2+(t47+t49+t51-t53)*t
     #24+(256.D0/3.D0*t69-256.D0/3.D0*t70)*t28+(256.D0/3.D0*t92-256.D0/3
     #.D0*s4l2)*s2l1+(-256.D0/3.D0*t96+256.D0/3.D0*s4l1)*s2l2+t65-t66)*t
     #19+((-256.D0/3.D0*s2l2-256.D0/3.D0*s4l2)*t104+((-512.D0/3.D0*s2l2+
     #512.D0/3.D0*s3l2+512.D0/3.D0*s4l2)*t28-t21+t26-256.D0/3.D0*t57)*t2
     #4+(-256.D0/3.D0*s2l2-256.D0/3.D0*s3l2)*t114+(-t21-256.D0/3.D0*t92+
     #t22)*t28+256.D0/3.D0*t121*s2l2-512.D0/3.D0*t124-512.D0/3.D0*t126)*
     #s1l1
      taux4 = s1+((256.D0/3.D0*s2l1+256.D0/3.D0*s4l1)*t104+((512.D0/3.D0
     #*s2l1-512.D0/3.D0*s3l1-512.D0/3.D0*s4l1)*t28+t33-t37+256.D0/3.D0*t
     #61)*t24+(256.D0/3.D0*s2l1+256.D0/3.D0*s3l1)*t114+(t33+256.D0/3.D0*
     #t96-t34)*t28-256.D0/3.D0*t121*s2l1+512.D0/3.D0*t147+512.D0/3.D0*t1
     #49)*s1l2+(-512.D0/3.D0*t69+512.D0/3.D0*t70+512.D0/3.D0*t50-512.D0/
     #3.D0*t52)*t24+(-512.D0/3.D0*t3+512.D0/3.D0*t4-512.D0/3.D0*t50+512.
     #D0/3.D0*t52)*t28+(512.D0/3.D0*t124+512.D0/3.D0*t126)*s2l1+(-512.D0
     #/3.D0*t147-512.D0/3.D0*t149)*s2l2

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux4

CCC Jet functions

      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet5 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet6 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff6)

CCC Sum

      rm2num2e2 = tFjet5*tcoff5+tFjet6*tcoff6
      RETURN
      END

      double precision function rm2num2e3(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin,s34spinc
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag,s34spinc
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      t1 = s24f(z,l1,l2,l3,l4)
      t2 = z*t1
      t5 = s13f(z,l1,l2,l3,l4)
      t6 = t5**2
      t8 = s23f(z,l1,l2,l3,l4)
      t9 = z*t8
      t11 = 64.D0/3.D0*t2
      t12 = 128.D0/3.D0*z
      t14 = s14f(z,l1,l2,l3,l4)
      t18 = t1**2
      t20 = 64.D0/3.D0*z*t18
      t22 = z*(-3+z)
      t27 = z*(-2+2*z)
      t28 = 64.D0/3.D0*t27
      t33 = t14**2
      t35 = t8**2
      t37 = 64.D0/3.D0*z*t35
      t41 = 128.D0/3.D0*t2
      t49 = z**2
      taux3 = (-64.D0/3.D0*t2-64.D0/3.D0*z)*t6+((-64.D0/3.D0*t9-t11-t12)
     #*t14+(-t11-t12)*t8-t20+64.D0/3.D0*t22*t1+t28)*t5+(-64.D0/3.D0*t9-6
     #4.D0/3.D0*z)*t33+(-t37+(-64.D0/3.D0*t2+64.D0/3.D0*t22)*t8-t41+t28)
     #*t14-t37+(-t41+t28)*t8-t20+64.D0/3.D0*t27*t1-64.D0/3.D0*z*(-2*z+t4
     #9+1)
      t1 = s1l1*s2l2
      t2 = s1l2*s2l1
      t3 = s2l1*s4l2
      t4 = s2l2*s4l1
      t7 = s13f(z,l1,l2,l3,l4)
      t8 = t7**2
      t19 = s14f(z,l1,l2,l3,l4)
      t21 = 2560.D0/3.D0*s2l2
      t22 = 1280.D0/3.D0*s4l2
      t24 = s23f(z,l1,l2,l3,l4)
      t26 = 1280.D0/3.D0*s3l2
      t27 = 512.D0/3.D0*s4l2
      t29 = s24f(z,l1,l2,l3,l4)
      t31 = -1+z
      t32 = 10*t31
      t34 = 256.D0/3.D0*t32*s2l2
      t37 = 2560.D0/3.D0*s2l1
      t38 = 1280.D0/3.D0*s4l1
      t41 = 1280.D0/3.D0*s3l1
      t42 = 512.D0/3.D0*s4l1
      t46 = 256.D0/3.D0*t32*s2l1
      t52 = 512.D0/3.D0*s3l2
      t55 = 512.D0/3.D0*s3l1
      t59 = 1280.D0/3.D0*s3l1*s4l2
      t61 = 1280.D0/3.D0*s3l2*s4l1
      t64 = 5*t31
      t66 = 256.D0/3.D0*t64*s4l2
      t70 = 256.D0/3.D0*t64*s4l1
      t75 = s2l1*s3l2
      t76 = s2l2*s3l1
      t79 = t19**2
      t103 = 256.D0/3.D0*t64*s3l2
      t107 = 256.D0/3.D0*t64*s3l1
      t114 = t24**2
      t124 = t29**2
      t129 = z**2
      t131 = -10*z+5*t129+5
      s1 = (1280.D0/3.D0*t1-1280.D0/3.D0*t2-1280.D0/3.D0*t3+1280.D0/3.D0
     #*t4)*t8+((2560.D0/3.D0*t1-2560.D0/3.D0*t2+(-1792.D0/3.D0*s3l2-1792
     #.D0/3.D0*s4l2)*s2l1+(1792.D0/3.D0*s3l1+1792.D0/3.D0*s4l1)*s2l2)*t1
     #9+((t21+t22)*t24+(t21+t26+t27)*t29-t34)*s1l1+((-t37-t38)*t24+(-t37
     #-t41-t42)*t29+t46)*s1l2+(-1280.D0/3.D0*t3+1280.D0/3.D0*t4)*t24+((-
     #t52-t22)*s2l1+(t55+t38)*s2l2-t59+t61)*t29+(-t52+t66)*s2l1+(t55-t70
     #)*s2l2)*t7+(1280.D0/3.D0*t1-1280.D0/3.D0*t2-1280.D0/3.D0*t75+1280.
     #D0/3.D0*t76)*t79
      taux4 = s1+(((t21+t52+t22)*t24+(t21+t26)*t29-t34)*s1l1+((-t37-t55-
     #t38)*t24+(-t37-t41)*t29+t46)*s1l2+((-t26-t27)*s2l1+(t41+t42)*s2l2+
     #t59-t61)*t24+(-1280.D0/3.D0*t75+1280.D0/3.D0*t76)*t29+(t103-t27)*s
     #2l1+(-t107+t42)*s2l2)*t19+((1280.D0/3.D0*s2l2+1280.D0/3.D0*s4l2)*t
     #114+((t21+1792.D0/3.D0*s3l2+1792.D0/3.D0*s4l2)*t29-t34+t52-t66)*t2
     #4+(1280.D0/3.D0*s2l2+1280.D0/3.D0*s3l2)*t124+(-t34-t103+t27)*t29+2
     #56.D0/3.D0*t131*s2l2)*s1l1+((-1280.D0/3.D0*s2l1-1280.D0/3.D0*s4l1)
     #*t114+((-t37-1792.D0/3.D0*s3l1-1792.D0/3.D0*s4l1)*t29+t46-t55+t70)
     #*t24+(-1280.D0/3.D0*s2l1-1280.D0/3.D0*s3l1)*t124+(t46+t107-t42)*t2
     #9-256.D0/3.D0*t131*s2l1)*s1l2

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux4

CCC Jet functions

      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet5 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet6 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff6)

CCC Sum

      rm2num2e3 = tFjet5*tcoff5+tFjet6*tcoff6
      RETURN
      END

      double precision function rm2num2e4(var,ndim,z,ux,l1,
     & l2,l3,l4,l5,l6,l7,l8,flag1,flag2,qflag,Tcoeff)
      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision lum,Pi,z3,NF
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,r,t,ux,Fjet
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,s34r2,
     & zs13f,zs14f,zs23f,zs24f
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf,Fjet,
     & zs13f,zs14f,zs23f,zs24f
      double precision s34r,s13r,s23r,s1vr,s3vr,s4vr,s2vr,
     & s134r,s234r,flag1,flag2
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      external fillspin,zfillspin,s34spinc
      double precision s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v
      double precision s34p,s13p,s23p,s1vp,s2vp,s14p,s24p,s3l1p,s4l1p,
     & s3l2p,s4l2p,qflag,s34spinc
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      if (flag1.lt.(0.0d0)) then
         s34r=s34rf(z,l1,l2,l3,l4)
         s13r=s13rf(z,l1,l2,l3,l4)
         s23r=s23rf(z,l1,l2,l3,l4)
         s1vr=s1vrf(z,l1,l2,l3,l4)
         s2vr=s1vrf(z,l1,l2,1-l3,l4)
         s234r=-(1+l1*(1-r*t)/r/(r+t))
         s1=r*(r+t)*(1+r*t)*(1+l1*(1-r*t)/r/(r+t))
         s2=(r+t)*(1+r*t)
         s3=(r-t)*(1-r*t)
         s134r=-(l1*(1+u)**2+(1-l1)*s2*r*(1-l1*s3/s2))/s1
         s34=l3*(1.0d0-l3)*(1.0d0-z)**2*s34r
         s13=(1.0d0-l3)*(1.0d0-z)*s13r
         s23=l3*(1.0d0-z)*s23r
         s1v=l3*(1.0d0-z)*s1vr
         s2v=(1.0d0-l3)*(1.0d0-z)*s2vr
         s14=-1.0d0-s13-s1v+z
         s24=s1v-s23-s34
         s3v=z-s13-s23-s34
         s4v=1.0d0+s13+s23
         s34p=l3*(1.0d0-l3)*(1.0d0-z)*s34r
         s13p=(1.0d0-l3)*s13r
         s23p=l3*s23r
         s1vp=l3*s1vr
         s2vp=(1.0d0-l3)*s2vr
         s14p=-1.0d0-s13p-s1vp
         s24p=s1vp-s23p-(1.0d0-z)*s34p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      else
         s13=s13f(z,l1,l2,l3,l4)
         s14=s14f(z,l1,l2,l3,l4)
         s23=s23f(z,l1,l2,l3,l4)
         s24=s24f(z,l1,l2,l3,l4)
         s34=-1.0d0+z-s13-s14-s23-s24
         s1v=-1.0d0+z-s13-s14
         s2v=-1.0d0+z-s23-s24
         s3v=1.0d0+s14+s24
         s4v=1.0d0+s13+s23
         s13p=zs13f(z,l1,l2,l3,l4)
         s14p=zs14f(z,l1,l2,l3,l4)
         s23p=zs23f(z,l1,l2,l3,l4)
         s24p=zs24f(z,l1,l2,l3,l4)
         s34p=-1.0d0-s13p-s14p-s23p-s24p
         s1vp=-1.0d0-s13p-s14p
         s2vp=-1.0d0-s23p-s24p
      call zfillspin(z,s13p,s14p,s23p,s24p,s34p,s1vp,s2vp,s3v,s4v,
     & l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1p,s4l1p,s1l2,s2l2,s3l2p,
     & s4l2p,c3l1,qflag)
      endif
      call fillspin(z,s13,s14,s23,s24,s34,s1v,s2v,s3v,s4v,l1,l2,l3,l4,
     & l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1,qflag)
      s3vr=z-(1-l3)*(1-z)*s13r-l3*(1-z)*s23r
     & -l3*(1-l3)*(1-z)**2*s34r
      s4vr=1.0d0+(1-l3)*(1-z)*s13r+l3*(1-z)*s23r

CCC Auxiliary functions

      taux3 = 0
      t1 = s1l1*s2l2
      t2 = s1l2*s2l1
      t3 = s2l1*s4l2
      t4 = s2l2*s4l1
      t7 = s13f(z,l1,l2,l3,l4)
      t8 = t7**2
      t19 = s14f(z,l1,l2,l3,l4)
      t21 = 1024.D0/3.D0*s2l2
      t22 = 512.D0/3.D0*s4l2
      t24 = s23f(z,l1,l2,l3,l4)
      t26 = 512.D0/3.D0*s3l2
      t28 = s24f(z,l1,l2,l3,l4)
      t30 = -1+z
      t31 = 2*t30
      t32 = t31*s2l2
      t33 = 512.D0/3.D0*t32
      t35 = ((t21+t22)*t24+(t21+t26)*t28-t33)*s1l1
      t36 = 1024.D0/3.D0*s2l1
      t37 = 512.D0/3.D0*s4l1
      t40 = 512.D0/3.D0*s3l1
      t43 = t31*s2l1
      t44 = 512.D0/3.D0*t43
      t46 = ((-t36-t37)*t24+(-t36-t40)*t28+t44)*s1l2
      t50 = s3l1*s4l2
      t51 = s3l2*s4l1
      t55 = t30*s4l2
      t58 = t30*s4l1
      t63 = s2l1*s3l2
      t64 = s2l2*s3l1
      t67 = t19**2
      t75 = t30*s3l2
      t78 = t30*s3l1
      t85 = t24**2
      t94 = t28**2
      t100 = z**2
      t101 = -2*z+t100+1
      s1 = (512.D0/3.D0*t1-512.D0/3.D0*t2-512.D0/3.D0*t3+512.D0/3.D0*t4)
     #*t8+((1024.D0/3.D0*t1-1024.D0/3.D0*t2+(-512.D0/3.D0*s3l2-512.D0/3.
     #D0*s4l2)*s2l1+(512.D0/3.D0*s3l1+512.D0/3.D0*s4l1)*s2l2)*t19+t35+t4
     #6+(-512.D0/3.D0*t3+512.D0/3.D0*t4)*t24+(-512.D0/3.D0*t3+512.D0/3.D
     #0*t4-512.D0/3.D0*t50+512.D0/3.D0*t51)*t28+512.D0/3.D0*t55*s2l1-512
     #.D0/3.D0*t58*s2l2)*t7+(512.D0/3.D0*t1-512.D0/3.D0*t2-512.D0/3.D0*t
     #63+512.D0/3.D0*t64)*t67
      taux4 = s1+(t35+t46+(-512.D0/3.D0*t63+512.D0/3.D0*t64+512.D0/3.D0*
     #t50-512.D0/3.D0*t51)*t24+(-512.D0/3.D0*t63+512.D0/3.D0*t64)*t28+51
     #2.D0/3.D0*t75*s2l1-512.D0/3.D0*t78*s2l2)*t19+((512.D0/3.D0*s2l2+51
     #2.D0/3.D0*s4l2)*t85+((t21+t26+t22)*t28-t33-512.D0/3.D0*t55)*t24+(5
     #12.D0/3.D0*s2l2+512.D0/3.D0*s3l2)*t94+(-512.D0/3.D0*t32-512.D0/3.D
     #0*t75)*t28+512.D0/3.D0*t101*s2l2)*s1l1+((-512.D0/3.D0*s2l1-512.D0/
     #3.D0*s4l1)*t85+((-t36-t40-t37)*t28+t44+512.D0/3.D0*t58)*t24+(-512.
     #D0/3.D0*s2l1-512.D0/3.D0*s3l1)*t94+(512.D0/3.D0*t43+512.D0/3.D0*t7
     #8)*t28-512.D0/3.D0*t101*s2l1)*s1l2

CCC Coefficients

      tcoff5 = taux3
      tcoff6 = taux4

CCC Jet functions

      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet5 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,3.D0,1.D0,1.D0,Tcoeff*tcoff5)
      t1 = s13f(z,l1,l2,l3,l4)
      t2 = s23f(z,l1,l2,l3,l4)
      t3 = s14f(z,l1,l2,l3,l4)
      t4 = s24f(z,l1,l2,l3,l4)
      tFjet6 = Fjet(var,ndim,z,ux,t1,t2,t3,t4,-1+z-t1-t3-t2-t4,-1+z-t1-t
     #3,-1+z-t2-t4,1+t3+t4,1+t1+t2,s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4
     #l2,6.D0,1.D0,1.D0,Tcoeff*tcoff6)

CCC Sum

      rm2num2e4 = tFjet5*tcoff5+tFjet6*tcoff6
      RETURN
      END

