CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgRV13e0int(var,ndim,Tcoeff)

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
      double precision qgRVnum9e0
      external qgRVnum9e0
      double precision qgRVnum9e1
      external qgRVnum9e1
      double precision qgRVnum9e2
      external qgRVnum9e2
      double precision qgRVnum9e3
      external qgRVnum9e3
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
      t7 = dlog(l1)
      t9 = t7**2
      t12 = dlog(z)
      t15 = t12*t2
      t21 = 0.3141592653589793D1**2
      t23 = t12**2
      t25 = 18*t7-3*t9-12*t3+12*t7*t12+24*t15-12*t7*t2-42-36*t12+36*t2+4
     #*t21-12*t23
      tncoff1 = -9.D0/512.D0*t3+t3*t2/256+t25/l1/2048+z3/128-t23*t12/256
     #+9.D0/256.D0*t15-45.D0/1024.D0-3.D0/256.D0*t3*t12-t2*t21/256+3.D0/
     #256.D0*t2*t23+t12*t21/256-21.D0/512.D0*t12+21.D0/512.D0*t2+3.D0/51
     #2.D0*t21-9.D0/512.D0*t23
      t1 = dlog(l1)
      t4 = dlog(1-z)
      t5 = t4**2
      t7 = dlog(z)
      t8 = t7**2
      t10 = 0.3141592653589793D1**2
      t13 = dlog(1-l1)
      t14 = t13**2
      t19 = t1**2
      t33 = -18*t1+12*t5+12*t8-4*t10+3*t14+36*t7-18*t13-36*t4+3*t19+42-2
     #4*t7*t4+12*t13*t4-12*t13*t7-12*t1*t7+12*t1*t4+6*t1*t13
      tncoff2 = t33/l1/2048
      t1 = dlog(z)
      t2 = t1**2
      t4 = 0.3141592653589793D1**2
      t7 = dlog(1-z)
      t8 = t7**2
      t14 = dlog(l1)
      tncoff3 = -3.D0/512.D0*t2+t4/512-3.D0/512.D0*t8-9.D0/512.D0*t1+9.D
     #0/512.D0*t7-21.D0/1024.D0+3.D0/256.D0*t1*t7+(6*t14-12*t1+12*t7-18)
     #/l1/2048
      t1 = dlog(z)
      t4 = dlog(1-l1)
      t7 = dlog(1-z)
      t9 = dlog(l1)
      tncoff4 = (12*t1-6*t4-12*t7-6*t9+18)/l1/2048
      t2 = dlog(1-z)
      t6 = dlog(z)
      tncoff5 = 3.D0/512.D0*t2-3.D0/1024.D0/l1-9.D0/1024.D0-3.D0/512.D0*
     #t6
      tncoff6 = 3.D0/1024.D0/l1
      tncoff7 = -3.D0/1024.D0

CCC Numerator functions

      tcnum1 = qgRVnum9e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum9e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum9e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum9e1(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum9e2(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgRVnum9e2(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgRVnum9e3(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff7)

CCC Sum

      qgRV13e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV13e1int(var,ndim,Tcoeff)

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
      double precision qgRVnum9e0
      external qgRVnum9e0
      double precision qgRVnum9e1
      external qgRVnum9e1
      double precision qgRVnum9e2
      external qgRVnum9e2
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
      t14 = dlog(l1)
      tncoff1 = -3.D0/512.D0*t3+9.D0/512.D0*t2+3.D0/256.D0*t6*t2+t9/512-
     #21.D0/1024.D0-9.D0/512.D0*t6-3.D0/512.D0*t12-3.D0/1024.D0*(-t14-2*
     #t2+3+2*t6)/l1
      t1 = dlog(l1)
      t3 = dlog(1-l1)
      t4 = dlog(z)
      t7 = dlog(1-z)
      tncoff2 = -3.D0/1024.D0*(t1+t3-3-2*t4+2*t7)/l1
      t2 = dlog(1-z)
      t6 = dlog(z)
      tncoff3 = 3.D0/512.D0*t2-3.D0/1024.D0/l1-9.D0/1024.D0-3.D0/512.D0*
     #t6
      tncoff4 = 3.D0/1024.D0/l1
      tncoff5 = -3.D0/1024.D0

CCC Numerator functions

      tcnum1 = qgRVnum9e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum9e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum9e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum9e1(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum9e2(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      qgRV13e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV13e2int(var,ndim,Tcoeff)

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
      double precision qgRVnum9e0
      external qgRVnum9e0
      double precision qgRVnum9e1
      external qgRVnum9e1
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
      t6 = dlog(z)
      tncoff1 = 3.D0/512.D0*t2-3.D0/1024.D0/l1-9.D0/1024.D0-3.D0/512.D0*
     #t6
      tncoff2 = 3.D0/1024.D0/l1
      tncoff3 = -3.D0/1024.D0

CCC Numerator functions

      tcnum1 = qgRVnum9e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum9e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum9e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff3)

CCC Sum

      qgRV13e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV13e3int(var,ndim,Tcoeff)

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
      double precision qgRVnum9e0
      external qgRVnum9e0
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

      tcnum1 = qgRVnum9e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      qgRV13e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

