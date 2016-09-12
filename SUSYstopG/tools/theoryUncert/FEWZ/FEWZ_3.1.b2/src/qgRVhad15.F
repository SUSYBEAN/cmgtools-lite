CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgRV15e0int(var,ndim,Tcoeff)

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
      double precision qgRVnum11e0
      external qgRVnum11e0
      double precision qgRVnum11e1
      external qgRVnum11e1
      double precision qgRVnum11e2
      external qgRVnum11e2
      double precision qgRVnum11e3
      external qgRVnum11e3
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
      t6 = 0.3141592653589793D1**2
      t8 = t1**2
      t11 = t1*t4
      t16 = t4**2
      t18 = 36*t1
      t19 = 36*t4
      t20 = 4*t6
      t21 = 24*t11
      t22 = 12*t8
      t23 = 12*t16
      t24 = dlog(l1)
      t26 = t24**2
      t32 = -t18+t19+t20+t21-42-t22-t23+36*t24-12*t26-24*t24*t4+24*t24*t
     #1
      t33 = 1/l1
      t43 = dlog(l3)
      t45 = 1/l3
      t54 = t43**2
      t56 = t18-t19-t20-t21+42+t22+t23-18*t43-12*t43*t1+12*t43*t4+3*t54
      tncoff1 = -21.D0/1024.D0*t1+21.D0/1024.D0*t4+3.D0/1024.D0*t6+3.D0/
     #512.D0*t4*t8+9.D0/512.D0*t11-45.D0/2048.D0+t1*t6/512-9.D0/1024.D0*
     #t8-9.D0/1024.D0*t16+t32*t33/2048-t4*t6/512-3.D0/512.D0*t16*t1+3.D0
     #/1024.D0*(-2*t1+2*t4-3+2*t24+t43)*t45*t33-t56*t45/4096+t16*t4/512+
     #z3/256-t8*t1/512
      t2 = dlog(1-z)
      t4 = dlog(z)
      t6 = dlog(l3)
      t7 = dlog(l1)
      t10 = 1/l3
      t24 = 0.3141592653589793D1**2
      t26 = t4**2
      t28 = t2**2
      t30 = t6**2
      t32 = 12*t4*t6-12*t6*t2+24*t4*t2+36*t2-36*t4+18*t6+4*t24-12*t26-12
     #*t28-3*t30-42
      tncoff2 = 3.D0/1024.D0*(-2*t2+2*t4-t6-2*t7+3)*t10/l1-t32*t10/4096
      t2 = dlog(1-l1)
      t4 = dlog(l1)
      t5 = dlog(z)
      t11 = dlog(1-z)
      t15 = t4**2
      t25 = t11**2
      t27 = t2**2
      t29 = 0.3141592653589793D1**2
      t31 = t5**2
      t33 = -18*t2-24*t4*t5+12*t4*t2+24*t4*t11-36*t4+12*t15+36*t5-36*t11
     #-12*t2*t5+12*t2*t11-24*t5*t11+12*t25+3*t27-4*t29+12*t31+42
      t34 = 1/l1
      t40 = dlog(l3)
      tncoff3 = t33*t34/2048+3.D0/1024.D0*(-t2-2*t4+2*t5-2*t11-t40+3)/l3
     #*t34
      t2 = dlog(1-l1)
      t4 = dlog(1-z)
      t6 = dlog(z)
      t8 = dlog(l3)
      t9 = dlog(l1)
      tncoff4 = 3.D0/1024.D0*(t2+2*t4-2*t6+t8+2*t9-3)/l3/l1
      t2 = dlog(1-z)
      t3 = 12*t2
      t4 = dlog(z)
      t5 = 12*t4
      t6 = dlog(l1)
      t9 = 1/l1
      t16 = 1/l3
      t19 = dlog(l3)
      t24 = t2**2
      t26 = 0.3141592653589793D1**2
      t28 = t4**2
      tncoff5 = -21.D0/2048.D0+(-18+t3-t5+12*t6)*t9/2048+9.D0/1024.D0*t2
     #-9.D0/1024.D0*t4+3.D0/512.D0*t4*t2-3.D0/1024.D0*t16*t9-(18-t3+t5-6
     #*t19)*t16/4096-3.D0/1024.D0*t24+t26/1024-3.D0/1024.D0*t28
      t1 = 1/l3
      t6 = dlog(1-z)
      t8 = dlog(z)
      t10 = dlog(l3)
      tncoff6 = 3.D0/1024.D0*t1/l1-(12*t6-12*t8+6*t10-18)*t1/4096
      t1 = dlog(l1)
      t4 = dlog(1-l1)
      t7 = dlog(1-z)
      t9 = dlog(z)
      t12 = 1/l1
      tncoff7 = (-12*t1-6*t4-12*t7+12*t9+18)*t12/2048+3.D0/1024.D0/l3*t1
     #2
      tncoff8 = -3.D0/1024.D0/l3/l1
      t6 = dlog(1-z)
      t8 = dlog(z)
      tncoff9 = -3.D0/1024.D0/l1-3.D0/2048.D0/l3+3.D0/1024.D0*t6-3.D0/10
     #24.D0*t8-9.D0/2048.D0
      tncoff10 = 3.D0/2048.D0/l3
      tncoff11 = 3.D0/1024.D0/l1
      tncoff12 = -3.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum11e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum11e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum11e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum11e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum11e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgRVnum11e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgRVnum11e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgRVnum11e1(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgRVnum11e2(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qgRVnum11e2(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qgRVnum11e2(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qgRVnum11e3(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      qgRV15e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV15e1int(var,ndim,Tcoeff)

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
      double precision qgRVnum11e0
      external qgRVnum11e0
      double precision qgRVnum11e1
      external qgRVnum11e1
      double precision qgRVnum11e2
      external qgRVnum11e2
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
      t14 = 2*t2
      t15 = 2*t6
      t16 = dlog(l3)
      t18 = 1/l3
      t21 = dlog(l1)
      t24 = 1/l1
      tncoff1 = -21.D0/2048.D0-3.D0/1024.D0*t3+9.D0/1024.D0*t2+3.D0/512.
     #D0*t6*t2+t9/1024-9.D0/1024.D0*t6-3.D0/1024.D0*t12+3.D0/2048.D0*(t1
     #4-3-t15+t16)*t18-3.D0/1024.D0*(-t14+3+t15-2*t21)*t24-3.D0/1024.D0*
     #t18*t24
      t2 = dlog(1-z)
      t4 = dlog(z)
      t6 = dlog(l3)
      t8 = 1/l3
      tncoff2 = 3.D0/2048.D0*(-2*t2+3+2*t4-t6)*t8+3.D0/1024.D0*t8/l1
      t2 = dlog(1-l1)
      t3 = dlog(z)
      t6 = dlog(1-z)
      t8 = dlog(l1)
      t11 = 1/l1
      tncoff3 = -3.D0/1024.D0*(t2-3-2*t3+2*t6+2*t8)*t11+3.D0/1024.D0/l3*
     #t11
      tncoff4 = -3.D0/1024.D0/l3/l1
      t6 = dlog(1-z)
      t8 = dlog(z)
      tncoff5 = -3.D0/1024.D0/l1-3.D0/2048.D0/l3+3.D0/1024.D0*t6-3.D0/10
     #24.D0*t8-9.D0/2048.D0
      tncoff6 = 3.D0/2048.D0/l3
      tncoff7 = 3.D0/1024.D0/l1
      tncoff8 = -3.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum11e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum11e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum11e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum11e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum11e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgRVnum11e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgRVnum11e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgRVnum11e2(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff8)

CCC Sum

      qgRV15e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV15e2int(var,ndim,Tcoeff)

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
      double precision qgRVnum11e0
      external qgRVnum11e0
      double precision qgRVnum11e1
      external qgRVnum11e1
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

      t6 = dlog(1-z)
      t8 = dlog(z)
      tncoff1 = -3.D0/1024.D0/l1-3.D0/2048.D0/l3+3.D0/1024.D0*t6-3.D0/10
     #24.D0*t8-9.D0/2048.D0
      tncoff2 = 3.D0/2048.D0/l3
      tncoff3 = 3.D0/1024.D0/l1
      tncoff4 = -3.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum11e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum11e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum11e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum11e1(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      qgRV15e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV15e3int(var,ndim,Tcoeff)

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
      double precision qgRVnum11e0
      external qgRVnum11e0
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

      tncoff1 = -3.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum11e0(var,ndim,z,ux,0.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      qgRV15e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

