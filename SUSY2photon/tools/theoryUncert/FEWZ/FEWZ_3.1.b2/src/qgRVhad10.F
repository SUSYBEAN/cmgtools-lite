CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgRV10e0int(var,ndim,Tcoeff)

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
      double precision qgRVnum7e0
      external qgRVnum7e0
      double precision qgRVnum7e1
      external qgRVnum7e1
      double precision qgRVnum7e2
      external qgRVnum7e2
      double precision qgRVnum7e3
      external qgRVnum7e3
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
      t3 = dlog(2.D0)
      t4 = t2*t3
      t5 = dlog(z)
      t8 = t3**2
      t11 = 0.3141592653589793D1**2
      t14 = dlog(l1)
      t15 = t14*t3
      t18 = t14**2
      t22 = t5**2
      t24 = t2**2
      t28 = t3*t5
      t31 = t2*t5
      t33 = t14*t5
      t41 = -84-24*t15*l1-24*t18+2*t11-24*t8-24*t22-54*t24+42*l1+72*t14+
     #24*t28*l1+72*t31-24*t33*l1-36*t31*l1-72*t5-48*t28+72*t4
      t52 = t14*t2
      t69 = 108*t2-72*t3+36*t3*l1-36*t14*l1-54*t2*l1+36*t5*l1+36*t52*l1-
     #36*t4*l1-72*t52+27*t24*l1+12*t22*l1-t11*l1+12*t8*l1+48*t33+12*t18*
     #l1+48*t15
      t90 = -3.D0/512.D0*t4*t5+t8*t5/512-t11/4096+3.D0/1024.D0*t8+(t41+t
     #69)/(-2+l1)/l1/6144+t8*t3/1536+3.D0/1024.D0*t22+t22*t5/1536-z3/768
     #+t11*t2/4096-9.D0/4096.D0*t24*t2+27.D0/4096.D0*t24-t11*t3/6144
      t109 = 15.D0/2048.D0+t22*t3/512-3.D0/1024.D0*t22*t2-9.D0/1024.D0*t
     #31-t11*t5/6144-3.D0/1024.D0*t2*t8+7.D0/1024.D0*t5+9.D0/2048.D0*t24
     #*t3+3.D0/512.D0*t28-9.D0/1024.D0*t4-21.D0/2048.D0*t2+7.D0/1024.D0*
     #t3+9.D0/2048.D0*t24*t5
      tncoff1 = t90+t109
      t1 = dlog(l1)
      t4 = dlog(1-z)
      t5 = dlog(2.D0)
      t8 = 2-l1
      t9 = dlog(t8)
      t12 = dlog(z)
      t27 = 84-72*t1-108*t4*t5-36*t5*t9+72*t5*t12+36*t4*t9-72*t4*t12-24*
     #t12*t9+72*t1*t4+24*t1*t9-72*t1*t5
      t30 = t9**2
      t32 = t4**2
      t36 = t12**2
      t38 = t5**2
      t40 = t1**2
      t44 = 0.3141592653589793D1**2
      t46 = -48*t1*t12+6*t30+54*t32+72*t12-36*t9+24*t36+54*t38+24*t40-10
     #8*t4+108*t5-2*t44
      tncoff2 = -(t27+t46)/t8/l1/6144
      t1 = dlog(z)
      t3 = dlog(2.D0)
      t5 = t3**2
      t8 = dlog(1-z)
      t9 = t8**2
      t11 = 0.3141592653589793D1**2
      t15 = t1**2
      t32 = dlog(l1)
      tncoff3 = 7.D0/2048.D0+3.D0/1024.D0*t1+3.D0/1024.D0*t3+t5/1024+9.D
     #0/4096.D0*t9-t11/12288-3.D0/1024.D0*t8*t3+t15/1024-9.D0/2048.D0*t8
     #+t3*t1/512-3.D0/1024.D0*t8*t1+(-36-24*t1-24*t3+36*t8+18*l1-18*t8*l
     #1+12*t1*l1+12*t3*l1+24*t32-12*t32*l1)/(-2+l1)/l1/6144
      t1 = dlog(2.D0)
      t4 = dlog(1-z)
      t6 = 2-l1
      t7 = dlog(t6)
      t9 = dlog(z)
      t11 = dlog(l1)
      tncoff4 = -(36+36*t1-36*t4-12*t7+24*t9-24*t11)/t6/l1/6144
      t1 = dlog(z)
      t12 = dlog(1-z)
      t14 = dlog(2.D0)
      tncoff5 = t1/1024+(6*l1-12)/(-2+l1)/l1/6144-3.D0/2048.D0*t12+t14/1
     #024+3.D0/2048.D0
      tncoff6 = 1/(-2+l1)/l1/512
      tncoff7 = 1.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum7e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum7e0(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum7e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum7e1(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum7e2(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgRVnum7e2(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgRVnum7e3(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff7)

CCC Sum

      qgRV10e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV10e1int(var,ndim,Tcoeff)

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
      double precision qgRVnum7e0
      external qgRVnum7e0
      double precision qgRVnum7e1
      external qgRVnum7e1
      double precision qgRVnum7e2
      external qgRVnum7e2
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

      t1 = dlog(2.D0)
      t2 = t1**2
      t5 = dlog(1-z)
      t9 = dlog(z)
      t12 = t5**2
      t17 = 0.3141592653589793D1**2
      t20 = t9**2
      t22 = dlog(l1)
      tncoff1 = t2/1024-3.D0/1024.D0*t1*t5+3.D0/1024.D0*t1+t1*t9/512+9.D
     #0/4096.D0*t12-9.D0/2048.D0*t5-3.D0/1024.D0*t5*t9+7.D0/2048.D0-t17/
     #12288+3.D0/1024.D0*t9+t20/1024-(-4*t22+2*t22*l1+4*t1-2*t1*l1-6*t5+
     #3*t5*l1+6-3*l1+4*t9-2*t9*l1)/(-2+l1)/l1/1024
      t1 = dlog(l1)
      t4 = dlog(1-z)
      t6 = 2-l1
      t7 = dlog(t6)
      t9 = dlog(2.D0)
      t11 = dlog(z)
      tncoff2 = (4*t1+6*t4+2*t7-6*t9-6-4*t11)/t6/l1/1024
      t1 = dlog(2.D0)
      t4 = dlog(1-z)
      t6 = dlog(z)
      t8 = 2-l1
      tncoff3 = t1/1024-3.D0/2048.D0*t4+3.D0/2048.D0+t6/1024+1/l1/1024
      tncoff4 = 1/(-2+l1)/l1/512
      tncoff5 = 1.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum7e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum7e0(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum7e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum7e1(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum7e2(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      qgRV10e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV10e2int(var,ndim,Tcoeff)

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
      double precision qgRVnum7e0
      external qgRVnum7e0
      double precision qgRVnum7e1
      external qgRVnum7e1
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
      t6 = dlog(1-z)
      t8 = dlog(2.D0)
      tncoff1 = t1/1024+1/l1/1024-3.D0/2048.D0*t6+t8/1024+3.D0/2048.D0
      tncoff2 = 1/(-2+l1)/l1/512
      tncoff3 = 1.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum7e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum7e0(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum7e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff3)

CCC Sum

      qgRV10e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV10e3int(var,ndim,Tcoeff)

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
      double precision qgRVnum7e0
      external qgRVnum7e0
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

      tncoff1 = 1.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum7e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      qgRV10e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

