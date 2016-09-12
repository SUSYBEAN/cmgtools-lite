CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgRV4e0int(var,ndim,Tcoeff)

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
      double precision qgRVnum3e0
      external qgRVnum3e0
      double precision qgRVnum3e1
      external qgRVnum3e1
      double precision qgRVnum3e2
      external qgRVnum3e2
      double precision qgRVnum3e3
      external qgRVnum3e3
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
      t5 = dlog(z)
      t6 = t2*t5
      t8 = dlog(2.D0)
      t9 = t2*t8
      t11 = t8*t5
      t13 = 0.3141592653589793D1**2
      t16 = dlog(l1)
      t17 = t16*t5
      t21 = t16**2
      t25 = t16*t2
      t33 = t16*t8
      t37 = -28-32*t3+16*t17+16*t3*l1-8*t21+4*t21*l1-32*t25+32*t6+32*t9-
     #16*t11-l1*t13+2*t13-24*t5+16*t33+48*t2-24*t8
      t56 = t8**2
      t61 = t5**2
      t68 = 14*l1-8*t17*l1+16*t25*l1+24*t16-8*t33*l1-16*t9*l1-16*t6*l1+8
     #*t11*l1-12*t16*l1+12*t8*l1+4*t56*l1+12*t5*l1+4*t61*l1-24*t2*l1-8*t
     #61-8*t56
      t89 = -3.D0/256.D0*t3+3.D0/256.D0*t6+3.D0/256.D0*t9-3.D0/512.D0*t1
     #1+3.D0/4096.D0*t13-(t37+t68)/(-2+l1)/l1/2048-t8*t61/512-t56*t5/512
     #+t2*t61/256+t2*t56/256-t3*t5/128-t3*t8/128-7.D0/1024.D0*t5
      t109 = -15.D0/2048.D0-t2*t13/1024+t13*t8/2048+t13*t5/2048+7.D0/512
     #.D0*t2-7.D0/1024.D0*t8+t9*t5/128-t56*t8/1536-t61*t5/1536+z3/768+t3
     #*t2/192-3.D0/1024.D0*t61-3.D0/1024.D0*t56
      tncoff1 = t89+t109
      t1 = dlog(2.D0)
      t2 = dlog(z)
      t6 = dlog(1-z)
      t9 = dlog(l1)
      t12 = 2-l1
      t13 = dlog(t12)
      t28 = t2**2
      t30 = 32*t1*t2-32*t6*t2-32*t9*t1+16*t9*t13-32*t1*t13-16*t2*t13-64*
     #t6*t1+32*t9*t6-16*t9*t2+32*t6*t13+8*t28
      t31 = 0.3141592653589793D1**2
      t33 = t6**2
      t35 = t13**2
      t41 = t1**2
      t43 = t9**2
      t46 = 28-2*t31+32*t33+8*t35-48*t6+24*t2-24*t13+48*t1+32*t41+8*t43-
     #24*t9
      tncoff2 = (t30+t46)/t12/l1/2048
      t1 = dlog(l1)
      t6 = dlog(2.D0)
      t9 = dlog(1-z)
      t11 = dlog(z)
      t26 = t11**2
      t28 = 0.3141592653589793D1**2
      t30 = t6**2
      t32 = t9**2
      tncoff3 = -(8*t1+6*l1-4*t1*l1-8*t6+16*t9-8*t11-8*t9*l1+4*t11*l1+4*
     #t6*l1-12)/(-2+l1)/l1/2048-t26/1024+t28/4096-t30/1024-t32/256-3.D0/
     #1024.D0*t6+3.D0/512.D0*t9-3.D0/1024.D0*t11+t9*t6/256+t9*t11/256-t6
     #*t11/512-7.D0/2048.D0
      t2 = dlog(1-z)
      t4 = dlog(z)
      t6 = 2-l1
      t7 = dlog(t6)
      t9 = dlog(2.D0)
      t11 = dlog(l1)
      tncoff4 = (-16*t2+8*t4-8*t7+16*t9-8*t11+12)/t6/l1/2048
      t1 = dlog(z)
      t3 = dlog(2.D0)
      t6 = dlog(1-z)
      tncoff5 = -t1/1024-t3/1024+t6/512-(2*l1-4)/(-2+l1)/l1/2048-3.D0/20
     #48.D0
      tncoff6 = -1/(-2+l1)/l1/512
      tncoff7 = -1.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum3e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum3e0(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum3e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum3e1(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum3e2(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgRVnum3e2(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgRVnum3e3(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff7)

CCC Sum

      qgRV4e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV4e1int(var,ndim,Tcoeff)

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
      double precision qgRVnum3e0
      external qgRVnum3e0
      double precision qgRVnum3e1
      external qgRVnum3e1
      double precision qgRVnum3e2
      external qgRVnum3e2
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
      t5 = dlog(2.D0)
      t9 = dlog(z)
      t12 = t5**2
      t17 = 0.3141592653589793D1**2
      t20 = t9**2
      t22 = dlog(l1)
      tncoff1 = -t3/256+t2*t5/256+3.D0/512.D0*t2+t9*t2/256-t12/1024-3.D0
     #/1024.D0*t5-t9*t5/512+t17/4096-7.D0/2048.D0-3.D0/1024.D0*t9-t20/10
     #24+(-4*t22+2*t22*l1-8*t2+4*t2*l1+4*t5-2*t5*l1+6-3*l1+4*t9-2*t9*l1)
     #/(-2+l1)/l1/1024
      t1 = dlog(l1)
      t3 = 2-l1
      t4 = dlog(t3)
      t7 = dlog(1-z)
      t9 = dlog(2.D0)
      t11 = dlog(z)
      tncoff2 = -(4*t1+4*t4+8*t7-8*t9-6-4*t11)/t3/l1/1024
      t2 = dlog(1-z)
      t4 = dlog(2.D0)
      t6 = dlog(z)
      t8 = 2-l1
      tncoff3 = t2/512-t4/1024-3.D0/2048.D0-t6/1024-1/l1/1024
      tncoff4 = -1/(-2+l1)/l1/512
      tncoff5 = -1.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum3e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum3e0(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum3e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum3e1(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum3e2(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      qgRV4e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV4e2int(var,ndim,Tcoeff)

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
      double precision qgRVnum3e0
      external qgRVnum3e0
      double precision qgRVnum3e1
      external qgRVnum3e1
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
      t4 = dlog(2.D0)
      t6 = dlog(z)
      tncoff1 = t2/512-t4/1024-3.D0/2048.D0-t6/1024-1/l1/1024
      tncoff2 = -1/(-2+l1)/l1/512
      tncoff3 = -1.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum3e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum3e0(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum3e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff3)

CCC Sum

      qgRV4e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV4e3int(var,ndim,Tcoeff)

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
      double precision qgRVnum3e0
      external qgRVnum3e0
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

      tncoff1 = -1.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum3e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      qgRV4e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

