CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qiqbic17e0int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac,e1temp
      double precision secsum,z3,NF,Li2
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      external fillspin,Li2
      double precision qiqbicnum17e0
      external qiqbicnum17e0
      double precision qiqbicnum17e1
      external qiqbicnum17e1
      double precision qiqbicnum17e2
      external qiqbicnum17e2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = var(7)
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = var(9)
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t1 = z-1
      t7 = dlog(z/(-l2*t1+z))
      t8 = t7**2
      t14 = 1/(-l2+l2*z-z)
      t16 = dlog(l2*t1*t14)
      t18 = 3.D0/256.D0*t1*t16
      t22 = t16**2
      t25 = 0.3141592653589793D1**2
      t29 = dlog(l1)
      t31 = 1/l1
      t34 = t1*t31
      t38 = dlog(-z*t14)
      tncoff1 = -3.D0/1024.D0*t1*t8+(t18-3.D0/512.D0*z+3.D0/512.D0)*t7-3
     #.D0/256.D0*t1*t22+t18+t1*(-12+t25)/2048+3.D0/512.D0*t29*t1*t31+3.D
     #0/256.D0*t34*t16-3.D0/512.D0*t34*t38-3.D0/512.D0*t34
      t1 = dlog(l1)
      t2 = z-1
      t4 = 1/l1
      t7 = t2*t4
      t11 = 1/(-l2+l2*z-z)
      t13 = dlog(l2*t2*t11)
      t17 = dlog(-z*t11)
      t21 = dlog(1-l1)
      tncoff2 = -3.D0/512.D0*t1*t2*t4-3.D0/256.D0*t7*t13+3.D0/512.D0*t7*
     #t17-3.D0/512.D0*t7*t21+3.D0/512.D0*t7
      t1 = 1-z
      t7 = dlog(z/(l2*t1+z))
      t9 = -t1
      t15 = dlog(l2*t9/(-l2+l2*z-z))
      tncoff3 = 3.D0/512.D0*t1*t7+3.D0/256.D0*t9*t15-3.D0/512.D0*z+3.D0/
     #512.D0-3.D0/512.D0*t9/l1
      tncoff4 = 3.D0/512.D0*(z-1)/l1
      tncoff5 = -3.D0/512.D0*z+3.D0/512.D0

CCC Numerator functions

      tcnum1 = qiqbicnum17e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqbicnum17e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqbicnum17e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qiqbicnum17e1(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qiqbicnum17e2(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      qiqbic17e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum
     #4*tncoff4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqbic17e1int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac,e1temp
      double precision secsum,z3,NF,Li2
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      external fillspin,Li2
      double precision qiqbicnum17e0
      external qiqbicnum17e0
      double precision qiqbicnum17e1
      external qiqbicnum17e1
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = var(7)
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = var(9)
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t1 = z-1
      t7 = dlog(z/(-l2*t1+z))
      t15 = dlog(l2*t1/(-l2+l2*z-z))
      tncoff1 = -3.D0/512.D0*t1*t7+3.D0/256.D0*t1*t15-3.D0/512.D0*z+3.D0
     #/512.D0-3.D0/512.D0*t1/l1
      tncoff2 = 3.D0/512.D0*(z-1)/l1
      tncoff3 = -3.D0/512.D0*z+3.D0/512.D0

CCC Numerator functions

      tcnum1 = qiqbicnum17e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qiqbicnum17e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qiqbicnum17e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff3)

CCC Sum

      qiqbic17e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qiqbic17e2int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac,e1temp
      double precision secsum,z3,NF,Li2
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      external fillspin,Li2
      double precision qiqbicnum17e0
      external qiqbicnum17e0
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = var(7)
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = var(9)
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      tncoff1 = -3.D0/512.D0*z+3.D0/512.D0

CCC Numerator functions

      tcnum1 = qiqbicnum17e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      qiqbic17e2int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

