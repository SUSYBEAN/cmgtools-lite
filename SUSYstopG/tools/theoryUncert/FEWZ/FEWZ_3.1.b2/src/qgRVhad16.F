CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgRV16e0int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision qgRVnum12e0
      external qgRVnum12e0
      double precision qgRVnum12e1
      external qgRVnum12e1
      double precision qgRVnum12e2
      external qgRVnum12e2
      double precision qgRVnum12e3
      external qgRVnum12e3
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l3 = jac*var(3)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t1 = dlog(z)
      t3 = dlog(1-z)
      t4 = t1*t3
      t6 = t3**2
      t9 = t1**2
      t13 = 0.3141592653589793D1**2
      t33 = 1-l1
      t34 = dlog(t33)
      t39 = t34**2
      t42 = -24*t4+12*t9+12*t6-36*t3-4*t13+36*t1+42+12*t34*t3-12*t34*t1+
     #3*t39-18*t34
      tncoff1 = 9.D0/256.D0*t4+t6*t3/256-9.D0/512.D0*t9-9.D0/512.D0*t6+2
     #1.D0/512.D0*t3+3.D0/512.D0*t13-21.D0/512.D0*t1-45.D0/1024.D0-3.D0/
     #256.D0*t6*t1-t3*t13/256+3.D0/256.D0*t3*t9+t1*t13/256+z3/128-t9*t1/
     #256-t42/t33/2048
      t1 = dlog(l1)
      t2 = t1**2
      t4 = 0.3141592653589793D1**2
      t6 = dlog(z)
      t7 = t6**2
      t11 = dlog(1-z)
      t14 = 1-l1
      t15 = dlog(t14)
      t16 = t15**2
      t19 = t11**2
      t33 = -3*t2+4*t4-12*t7+18*t1+36*t11-36*t6-42-3*t16+18*t15-12*t19-1
     #2*t15*t11+24*t6*t11-12*t1*t11-6*t1*t15+12*t15*t6+12*t6*t1
      tncoff2 = -t33/t14/2048
      t2 = dlog(1-z)
      t3 = t2**2
      t5 = 0.3141592653589793D1**2
      t7 = dlog(z)
      t8 = t7**2
      t12 = 1-l1
      t13 = dlog(t12)
      tncoff3 = -3.D0/512.D0*t3+t5/512-3.D0/512.D0*t8-(-12*t2+12*t7-6*t1
     #3+18)/t12/2048+3.D0/256.D0*t7*t2+9.D0/512.D0*t2-9.D0/512.D0*t7-21.
     #D0/1024.D0
      t1 = 1-l1
      t2 = dlog(t1)
      t4 = dlog(l1)
      t6 = dlog(z)
      t9 = dlog(1-z)
      tncoff4 = -(6*t2+6*t4-12*t6+12*t9-18)/t1/2048
      t1 = dlog(z)
      t4 = dlog(1-z)
      tncoff5 = -3.D0/512.D0*t1+3.D0/512.D0*t4+3.D0/1024.D0/(-1+l1)-9.D0
     #/1024.D0
      tncoff6 = -3.D0/1024.D0/(-1+l1)
      tncoff7 = -3.D0/1024.D0

CCC Numerator functions

      tcnum1 = qgRVnum12e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum12e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum12e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum12e1(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum12e2(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgRVnum12e2(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgRVnum12e3(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff7)

CCC Sum

      qgRV16e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV16e1int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision qgRVnum12e0
      external qgRVnum12e0
      double precision qgRVnum12e1
      external qgRVnum12e1
      double precision qgRVnum12e2
      external qgRVnum12e2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l3 = jac*var(3)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t2 = dlog(1-z)
      t3 = t2**2
      t6 = dlog(z)
      t9 = 0.3141592653589793D1**2
      t12 = t6**2
      t14 = 1-l1
      t15 = dlog(t14)
      tncoff1 = -3.D0/512.D0*t3+9.D0/512.D0*t2+3.D0/256.D0*t6*t2-21.D0/1
     #024.D0+t9/512-9.D0/512.D0*t6-3.D0/512.D0*t12+3.D0/1024.D0*(t15+2*t
     #2-3-2*t6)/t14
      t1 = 1-l1
      t2 = dlog(t1)
      t4 = dlog(1-z)
      t6 = dlog(z)
      t8 = dlog(l1)
      tncoff2 = 3.D0/1024.D0*(-t2-2*t4+3+2*t6-t8)/t1
      t1 = dlog(z)
      t4 = dlog(1-z)
      tncoff3 = -3.D0/512.D0*t1+3.D0/512.D0*t4+3.D0/1024.D0/(-1+l1)-9.D0
     #/1024.D0
      tncoff4 = -3.D0/1024.D0/(-1+l1)
      tncoff5 = -3.D0/1024.D0

CCC Numerator functions

      tcnum1 = qgRVnum12e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum12e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum12e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum12e1(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum12e2(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      qgRV16e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV16e2int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision qgRVnum12e0
      external qgRVnum12e0
      double precision qgRVnum12e1
      external qgRVnum12e1
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l3 = jac*var(3)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t1 = dlog(z)
      t4 = dlog(1-z)
      tncoff1 = -3.D0/512.D0*t1+3.D0/512.D0*t4+3.D0/1024.D0/(-1+l1)-9.D0
     #/1024.D0
      tncoff2 = -3.D0/1024.D0/(-1+l1)
      tncoff3 = -3.D0/1024.D0

CCC Numerator functions

      tcnum1 = qgRVnum12e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum12e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum12e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff3)

CCC Sum

      qgRV16e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV16e3int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision qgRVnum12e0
      external qgRVnum12e0
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l3 = jac*var(3)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      tncoff1 = -3.D0/1024.D0

CCC Numerator functions

      tcnum1 = qgRVnum12e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      qgRV16e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

