CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function RV7e0int(var,ndim,Tcoeff)

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
      t12 = 1-z
      t13 = dlog(t12)
      t18 = t13**2
      t23 = -1/t12
      t26 = dlog(l1)
      t27 = t26*l1
      t34 = t26**2
      t36 = t1*l1
      t38 = t26*t1
      t47 = -12*t27-6*t2-24*t1+2*t6+24*t26+24*l1-6*t34+12*t36+12*t38+3*t
     #34*l1+3*t2*l1-t6*l1-48-6*t38*l1
      t49 = 1/(l1-2)
      t51 = 1/l1
      tncoff1 = -t2*t1/6912-t2/1152+t6*t1/6912-t1/288-1.D0/144.D0+z3/432
     #+t6/3456+(-18*t13*t1+3*t2+12*t1-t6+27*t18-36*t13+24)*t23/2304-t47*
     #t49*t51/6912+(-t27-2*t1+6*t13+2*t26+2*l1+t36-3*t13*l1-4)*t49*t51*t
     #23/384
      t1 = dlog(l1)
      t3 = dlog(2.D0)
      t5 = 2-l1
      t6 = dlog(t5)
      t14 = t1**2
      t16 = 0.3141592653589793D1**2
      t18 = t3**2
      t20 = t6**2
      t22 = -24*t1+72*t3-48*t6+24*t6*t1-36*t1*t3-72*t3*t6+6*t14-2*t16+54
     #*t18+24*t20+48
      t24 = -1/t5
      t26 = 1/l1
      t32 = 1-z
      t33 = dlog(t32)
      tncoff2 = -t22*t24*t26/6912-(-2*t1+6*t3-4*t6-6*t33+4)*t24*t26/t32/
     #384
      t1 = 0.3141592653589793D1**2
      t2 = dlog(z)
      t3 = t2**2
      t5 = 1-z
      t6 = dlog(t5)
      t7 = dlog(2.D0)
      t17 = t7**2
      t19 = t6**2
      t21 = t1-12*t3+18*t6*t7+36*t6*t2-12*t7*t2-12*t7-24*t2+36*t6-3*t17-
     #27*t19-24
      t23 = -1/t5
      t27 = dlog(l1)
      tncoff3 = t21*t23/2304+(-2*l1-2*t27-2*t2*l1+t27*l1+3*t6*l1-t7*l1+2
     #*t7+4*t2-6*t6+4)/(l1-2)/l1*t23/384
      t1 = dlog(z)
      t3 = 2-l1
      t4 = dlog(t3)
      t6 = 1-z
      t7 = dlog(t6)
      t9 = dlog(l1)
      t11 = dlog(2.D0)
      tncoff4 = (-4*t1+4*t4+6*t7+2*t9-6*t11-4)/t3/l1/t6/384
      t1 = dlog(2.D0)
      t2 = t1**2
      t5 = 0.3141592653589793D1**2
      t8 = 1-z
      t9 = dlog(t8)
      t13 = -1/t8
      t18 = dlog(l1)
      t28 = 1/l1
      tncoff5 = -t2/2304-t1/576-1.D0/288.D0+t5/6912+(6*t1+12-18*t9)*t13/
     #2304-(-12*t1-24+12*l1+12*t18-6*t18*l1+6*t1*l1)/(l1-2)*t28/6912+t28
     #*t13/384
      t1 = dlog(l1)
      t3 = dlog(2.D0)
      t5 = 2-l1
      t6 = dlog(t5)
      t10 = -1/t5
      t12 = 1/l1
      tncoff6 = -(24-12*t1+36*t3-24*t6)*t10*t12/6912+t10*t12/(-1+z)/192
      t1 = 1-z
      t2 = dlog(t1)
      t4 = dlog(2.D0)
      t6 = dlog(z)
      t10 = -1/t1
      t13 = 2-l1
      tncoff7 = (18*t2-6*t4-12*t6-12)*t10/2304-1/l1*t10/384
      tncoff8 = -1/(l1-2)/l1/(-1+z)/192
      t1 = dlog(2.D0)
      tncoff9 = -t1/1152-1.D0/576.D0+1/(-1+z)/384-(-12+6*l1)/(l1-2)/l1/6
     #912
      tncoff10 = -1/(l1-2)/l1/576
      tncoff11 = -1/(-1+z)/384
      tncoff12 = -1.D0/1152.D0

CCC Numerator functions

      tcnum1 = RVnum4e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum4e0(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum4e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum4e0(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum4e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum4e1(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum4e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum4e1(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = RVnum4e2(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = RVnum4e2(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = RVnum4e2(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = RVnum4e3(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      RV7e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV7e1int(var,ndim,Tcoeff)

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
      t7 = 1-z
      t8 = dlog(t7)
      t12 = -1/t7
      t15 = dlog(l1)
      t25 = 1/l1
      tncoff1 = -t2/2304-t1/576-1.D0/288.D0+t5/6912+(t1+2-3*t8)*t12/384-
     #(-t15*l1+2*t15+t1*l1-2*t1+2*l1-4)/(l1-2)*t25/1152+t25*t12/384
      t1 = dlog(l1)
      t3 = dlog(2.D0)
      t5 = 2-l1
      t6 = dlog(t5)
      t10 = -1/t5
      t12 = 1/l1
      tncoff2 = -(-2*t1+6*t3-4*t6+4)*t10*t12/1152+t10*t12/(-1+z)/192
      t1 = dlog(2.D0)
      t2 = dlog(z)
      t4 = 1-z
      t5 = dlog(t4)
      t9 = -1/t4
      t11 = 2-l1
      tncoff3 = (-t1-2-2*t2+3*t5)*t9/384-1/l1*t9/384
      tncoff4 = -1/(l1-2)/l1/(-1+z)/192
      t1 = dlog(2.D0)
      tncoff5 = -t1/1152-1.D0/576.D0+1/(-1+z)/384-1/l1/1152
      tncoff6 = -1/(l1-2)/l1/576
      tncoff7 = -1/(-1+z)/384
      tncoff8 = -1.D0/1152.D0

CCC Numerator functions

      tcnum1 = RVnum4e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum4e0(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum4e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum4e0(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum4e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum4e1(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum4e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum4e2(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff8)

CCC Sum

      RV7e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV7e2int(var,ndim,Tcoeff)

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
      tncoff1 = -t1/1152-1.D0/576.D0+1/(-1+z)/384-1/l1/1152
      tncoff2 = -1/(l1-2)/l1/576
      tncoff3 = -1/(-1+z)/384
      tncoff4 = -1.D0/1152.D0

CCC Numerator functions

      tcnum1 = RVnum4e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum4e0(var,ndim,1.D0,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum4e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum4e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      RV7e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV7e3int(var,ndim,Tcoeff)

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

      tncoff1 = -1.D0/1152.D0

CCC Numerator functions

      tcnum1 = RVnum4e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      RV7e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV7e4int(var,ndim,Tcoeff)

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

      RV7e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

