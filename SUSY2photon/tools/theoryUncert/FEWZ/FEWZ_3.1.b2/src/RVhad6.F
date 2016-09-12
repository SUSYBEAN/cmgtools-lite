CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function RV6e0int(var,ndim,Tcoeff)

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
      double precision RVnum4e0
      external RVnum4e0
      double precision RVnum4e1
      external RVnum4e1
      double precision RVnum4e2
      external RVnum4e2
      double precision RVnum4e3
      external RVnum4e3
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
      t6 = 0.3141592653589793D1**2
      t14 = 1-z
      t15 = dlog(t14)
      t16 = t15**2
      t23 = -1/t14
      t30 = dlog(l1)
      t31 = t30**2
      t34 = t30*t1
      t36 = t30*l1
      t38 = t1*l1
      t47 = -24*t2-48*t1+2*t6-t6*l1+12*t31*l1+48*t34-24*t36+24*t38+12*t2
     #*l1-24*t31+48*t30+24*l1-24*t34*l1-48
      t49 = 1/(l1-2)
      t51 = 1/l1
      tncoff1 = -1.D0/288.D0-t2*t1/1728-t2/576+t6*t1/6912-t1/288+z3/864+
     #t6/6912+(12*t2+24*t1-t6+27*t16-36*t15-36*t15*t1+24)*t23/4608-t47*t
     #49*t51/6912+(-4*t1+6*t15-2*t36+2*t38-3*t15*l1+4*t30+2*l1-4)*t49*t5
     #1*t23/384
      t1 = dlog(2.D0)
      t2 = 2-l1
      t3 = dlog(t2)
      t6 = dlog(l1)
      t11 = 0.3141592653589793D1**2
      t13 = t1**2
      t15 = t3**2
      t17 = t6**2
      t22 = -36*t1*t3-72*t6*t1+24*t6*t3-2*t11+54*t13+6*t15+24*t17-24*t3-
     #48*t6+72*t1+48
      t24 = -1/t2
      t26 = 1/l1
      t32 = 1-z
      t33 = dlog(t32)
      tncoff2 = -t22*t24*t26/6912-(-2*t3-4*t6+6*t1-6*t33+4)*t24*t26/t32/
     #384
      t1 = 1-z
      t2 = dlog(t1)
      t3 = dlog(z)
      t6 = dlog(2.D0)
      t11 = t3**2
      t13 = t6**2
      t15 = 0.3141592653589793D1**2
      t16 = t2**2
      t21 = 36*t2*t3+36*t2*t6-24*t6*t3-12*t11-12*t13+t15-27*t16-24*t6-24
     #*t3+36*t2-24
      t23 = -1/t1
      t31 = dlog(l1)
      tncoff3 = t21*t23/4608+(4*t6+4*t3-6*t2+3*t2*l1-4*t31-2*l1+2*t31*l1
     #-2*t3*l1-2*t6*l1+4)/(l1-2)/l1*t23/384
      t1 = 1-z
      t2 = dlog(t1)
      t4 = dlog(l1)
      t6 = dlog(2.D0)
      t8 = 2-l1
      t9 = dlog(t8)
      t11 = dlog(z)
      tncoff4 = (-4+6*t2+4*t4-6*t6+2*t9-4*t11)/t8/l1/t1/384
      t1 = dlog(2.D0)
      t2 = t1**2
      t5 = 0.3141592653589793D1**2
      t8 = 1-z
      t9 = dlog(t8)
      t13 = -1/t8
      t19 = dlog(l1)
      t28 = 1/l1
      tncoff5 = -1.D0/576.D0-t2/1152-t1/576+t5/13824+(12+12*t1-18*t9)*t1
     #3/4608-(-24-24*t1+12*t1*l1+24*t19+12*l1-12*t19*l1)/(l1-2)*t28/6912
     #+t28*t13/384
      t1 = dlog(2.D0)
      t3 = 2-l1
      t4 = dlog(t3)
      t6 = dlog(l1)
      t10 = -1/t3
      t12 = 1/l1
      tncoff6 = -(36*t1-12*t4-24*t6+24)*t10*t12/6912+t10*t12/(-1+z)/192
      t1 = dlog(2.D0)
      t3 = dlog(z)
      t5 = 1-z
      t6 = dlog(t5)
      t10 = -1/t5
      t13 = 2-l1
      tncoff7 = (-12*t1-12*t3+18*t6-12)*t10/4608-1/l1*t10/384
      tncoff8 = -1/(l1-2)/l1/(-1+z)/192
      t1 = dlog(2.D0)
      tncoff9 = -1.D0/1152.D0-t1/1152+1/(-1+z)/768-(-12+6*l1)/(l1-2)/l1/
     #6912
      tncoff10 = -1/(l1-2)/l1/576
      tncoff11 = -1/(-1+z)/768
      tncoff12 = -1.D0/2304.D0

CCC Numerator functions

      tcnum1 = RVnum4e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum4e0(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum4e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum4e0(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum4e1(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum4e1(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum4e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum4e1(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = RVnum4e2(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = RVnum4e2(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = RVnum4e2(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = RVnum4e3(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      RV6e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV6e1int(var,ndim,Tcoeff)

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
      double precision RVnum4e0
      external RVnum4e0
      double precision RVnum4e1
      external RVnum4e1
      double precision RVnum4e2
      external RVnum4e2
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
      t5 = 0.3141592653589793D1**2
      t8 = 1-z
      t9 = dlog(t8)
      t13 = -1/t8
      t16 = dlog(l1)
      t28 = 1/l1
      tncoff1 = -1.D0/576.D0-t2/1152-t1/576+t5/13824+(2+2*t1-3*t9)*t13/7
     #68-(-2*t16*l1+4*t16+2*t1*l1-4*t1+2*l1-4)/(l1-2)*t28/1152+t28*t13/3
     #84
      t1 = dlog(l1)
      t3 = dlog(2.D0)
      t5 = 2-l1
      t6 = dlog(t5)
      t10 = -1/t5
      t12 = 1/l1
      tncoff2 = -(-4*t1+6*t3-2*t6+4)*t10*t12/1152+t10*t12/(-1+z)/192
      t1 = dlog(2.D0)
      t3 = dlog(z)
      t5 = 1-z
      t6 = dlog(t5)
      t10 = -1/t5
      t13 = 2-l1
      tncoff3 = (-2*t1-2*t3-2+3*t6)*t10/768-1/l1*t10/384
      tncoff4 = -1/(l1-2)/l1/(-1+z)/192
      t1 = dlog(2.D0)
      tncoff5 = -1.D0/1152.D0-t1/1152+1/(-1+z)/768-1/l1/1152
      tncoff6 = -1/(l1-2)/l1/576
      tncoff7 = -1/(-1+z)/768
      tncoff8 = -1.D0/2304.D0

CCC Numerator functions

      tcnum1 = RVnum4e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum4e0(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum4e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum4e0(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum4e1(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum4e1(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum4e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum4e2(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff8)

CCC Sum

      RV6e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV6e2int(var,ndim,Tcoeff)

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
      double precision RVnum4e0
      external RVnum4e0
      double precision RVnum4e1
      external RVnum4e1
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
      tncoff1 = -1.D0/1152.D0-t1/1152+1/(-1+z)/768-1/l1/1152
      tncoff2 = -1/(l1-2)/l1/576
      tncoff3 = -1/(-1+z)/768
      tncoff4 = -1.D0/2304.D0

CCC Numerator functions

      tcnum1 = RVnum4e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum4e0(var,ndim,1.D0,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum4e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum4e1(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      RV6e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV6e3int(var,ndim,Tcoeff)

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
      double precision RVnum4e0
      external RVnum4e0
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

      tncoff1 = -1.D0/2304.D0

CCC Numerator functions

      tcnum1 = RVnum4e0(var,ndim,1.D0,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      RV6e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV6e4int(var,ndim,Tcoeff)

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

      RV6e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

