CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function RV9e0int(var,ndim,Tcoeff)

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
      double precision RVnum6e0
      external RVnum6e0
      double precision RVnum6e1
      external RVnum6e1
      double precision RVnum6e2
      external RVnum6e2
      double precision RVnum6e3
      external RVnum6e3
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

      t1 = 0.3141592653589793D1**2
      t5 = 1-z
      t6 = dlog(t5)
      t7 = t6**2
      t12 = -1/t5
      t15 = dlog(l1)
      t18 = t15**2
      t21 = 1/l1
      tncoff1 = 1.D0/32.D0-t1/192-z3/96+(2*t1-12-6*t7+12*t6)*t12/384-(12
     #*t15+4*t1-24-3*t18)*t21/1536+(2*t6-2+t15)*t21*t12/128
      t1 = dlog(l1)
      t3 = dlog(1-l1)
      t7 = t3**2
      t10 = 0.3141592653589793D1**2
      t12 = t1**2
      t15 = 1/l1
      t18 = 1-z
      t19 = dlog(t18)
      tncoff2 = -(6*t1*t3-12*t1+3*t7-12*t3-4*t10+24+3*t12)*t15/1536-(-2*
     #t19-t3+2-t1)*t15/t18/128
      t1 = 0.3141592653589793D1**2
      t3 = dlog(z)
      t5 = t3**2
      t7 = 1-z
      t8 = dlog(t7)
      t9 = t8**2
      t16 = -1/t7
      t20 = dlog(l1)
      tncoff3 = (-2*t1+12*t3+12+6*t5+6*t9-12*t8-12*t8*t3)*t16/384+(-2*t8
     #+2-t20+2*t3)/l1*t16/128
      t1 = 1-z
      t2 = dlog(t1)
      t5 = dlog(1-l1)
      t6 = dlog(z)
      t8 = dlog(l1)
      tncoff4 = -(2*t2+t5-2*t6+t8-2)/l1/t1/128
      t1 = 0.3141592653589793D1**2
      t3 = 1-z
      t4 = dlog(t3)
      t8 = -1/t3
      t11 = dlog(l1)
      t14 = 1/l1
      tncoff5 = -t1/384+1.D0/64.D0+(-6+6*t4)*t8/384-(6*t11-12)*t14/1536-
     #t14*t8/128
      t1 = dlog(l1)
      t4 = dlog(1-l1)
      t7 = 1/l1
      tncoff6 = -(-6*t1+12-6*t4)*t7/1536+t7/(-1+z)/128
      t1 = dlog(z)
      t2 = 1-z
      t3 = dlog(t2)
      t7 = -1/t2
      tncoff7 = (6*t1+6-6*t3)*t7/384+1/l1*t7/128
      tncoff8 = -1/l1/(-1+z)/128
      tncoff9 = 1.D0/128.D0-1/(-1+z)/128+1/l1/256
      tncoff10 = -1/l1/256
      tncoff11 = 1/(-1+z)/128
      tncoff12 = 1.D0/256.D0

CCC Numerator functions

      tcnum1 = RVnum6e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum6e0(var,ndim,1.D0,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum6e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum6e0(var,ndim,z,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum6e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum6e1(var,ndim,1.D0,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum6e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum6e1(var,ndim,z,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = RVnum6e2(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = RVnum6e2(var,ndim,1.D0,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = RVnum6e2(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = RVnum6e3(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      RV9e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV9e1int(var,ndim,Tcoeff)

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
      double precision RVnum6e0
      external RVnum6e0
      double precision RVnum6e1
      external RVnum6e1
      double precision RVnum6e2
      external RVnum6e2
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

      t1 = 0.3141592653589793D1**2
      t3 = 1-z
      t4 = dlog(t3)
      t8 = -1/t3
      t11 = dlog(l1)
      t13 = 1/l1
      tncoff1 = -t1/384+1.D0/64.D0+(-2+2*t4)*t8/128+(-t11+2)*t13/256-t13
     #*t8/128
      t1 = dlog(l1)
      t3 = dlog(1-l1)
      t5 = 1/l1
      tncoff2 = (t1+t3-2)*t5/256+t5/(-1+z)/128
      t1 = dlog(z)
      t2 = 1-z
      t3 = dlog(t2)
      t7 = -1/t2
      tncoff3 = (2+2*t1-2*t3)*t7/128+1/l1*t7/128
      tncoff4 = -1/l1/(-1+z)/128
      tncoff5 = 1.D0/128.D0-1/(-1+z)/128+1/l1/256
      tncoff6 = -1/l1/256
      tncoff7 = 1/(-1+z)/128
      tncoff8 = 1.D0/256.D0

CCC Numerator functions

      tcnum1 = RVnum6e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum6e0(var,ndim,1.D0,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum6e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum6e0(var,ndim,z,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum6e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum6e1(var,ndim,1.D0,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum6e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum6e2(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff8)

CCC Sum

      RV9e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV9e2int(var,ndim,Tcoeff)

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
      double precision RVnum6e0
      external RVnum6e0
      double precision RVnum6e1
      external RVnum6e1
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

      tncoff1 = 1.D0/128.D0-1/(-1+z)/128+1/l1/256
      tncoff2 = -1/l1/256
      tncoff3 = 1/(-1+z)/128
      tncoff4 = 1.D0/256.D0

CCC Numerator functions

      tcnum1 = RVnum6e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum6e0(var,ndim,1.D0,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum6e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum6e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      RV9e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV9e3int(var,ndim,Tcoeff)

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
      double precision RVnum6e0
      external RVnum6e0
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

      tncoff1 = 1.D0/256.D0

CCC Numerator functions

      tcnum1 = RVnum6e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      RV9e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV9e4int(var,ndim,Tcoeff)

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


CCC Numerator functions


CCC Sum

      RV9e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

