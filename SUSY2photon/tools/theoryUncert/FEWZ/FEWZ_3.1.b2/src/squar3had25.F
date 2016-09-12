CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function sq35e0int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
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
      double precision sq3num3e0
      external sq3num3e0
      double precision sq3num3e1
      external sq3num3e1
      double precision sq3num3e2
      external sq3num3e2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l4 = jac*l4+buf
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
      t3 = dlog(2.D0)
      t4 = t3**2
      t6 = dlog(l4)
      t10 = dlog(1-z)
      t11 = dlog(z)
      t23 = dlog(1-l4)
      t32 = f2(0.D0,z)
      t33 = dlog(t32)
      t38 = t1/6144-t4/256-t6*t3/512+t10*t11/256-t6*t10/512+t6*t11/1024-
     #t10*t3/128+t11*t3/256+t23*t11/1024-t23*t10/512-t23*t6/2048-t23*t3/
     #512+t23*t33/2048+t6*t33/2048
      t43 = dlog(l1)
      t58 = t11**2
      t60 = t33**2
      t62 = t10**2
      t64 = t23**2
      t66 = t6**2
      t68 = t43**2
      t70 = t10*t33/512+t33*t3/512-t23*t43/2048-t43*t10/512+t43*t11/1024
     #-t43*t3/512+t43*t33/2048-t6*t43/2048-t33*t11/1024-t58/1024-t60/409
     #6-t62/256-t64/4096-t66/4096-t68/4096
      t72 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t73 = 1/t72
      t75 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t77 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t79 = 1/t75/t77
      t80 = f1(0.D0,z,0.D0)
      t81 = t79*t80
      t83 = 4*t10
      t84 = 4*t3
      t85 = 2*t11
      t86 = dlog(l3)
      t89 = 1/l3
      t94 = dlog(l2)
      t97 = 1/l2
      tncoff1 = (t38+t70)*t73*t81+(t6+t23+t43+t83-t33+t84-t85+2*t86)*t89
     #*t73*t81/1024+(t23+t6+2*t94+t43+t83-t33+t84-t85)*t97*t73*t81/1024-
     #t79*t73*t97*t89*t80/512
      t1 = dlog(l4)
      t3 = dlog(1-l4)
      t4 = dlog(l1)
      t6 = dlog(1-l3)
      t9 = dlog(1-z)
      t11 = f2(l3,z)
      t12 = dlog(t11)
      t13 = dlog(2.D0)
      t15 = dlog(z)
      t17 = dlog(l3)
      t20 = 1/l3
      t22 = Dsquar321(z,l1,0.D0,l3,l4)
      t23 = 1/t22
      t25 = s23r2(z,0.D0,0.D0,l3,l4)
      t27 = s1vr2(z,0.D0,0.D0,l3,l4)
      t29 = 1/t25/t27
      t30 = f1(l3,z,0.D0)
      tncoff2 = (-t1-t3-t4-2*t6-4*t9+t12-4*t13+2*t15-2*t17)*t20*t23*t29*
     #t30/1024+t29*t23/l2*t20*t30/512
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = dlog(l2)
      t6 = dlog(l1)
      t8 = dlog(1-z)
      t11 = dlog(1-l2)
      t12 = l1*l2
      t14 = dlog(1-t12)
      t15 = f3(0.D0,z)
      t19 = dlog(1-t15*l1*l2)
      t20 = f2(0.D0,z)
      t21 = dlog(t20)
      t22 = dlog(2.D0)
      t24 = dlog(z)
      t26 = t2+t3+2*t4+t6+4*t8+t11+t14+t19-t21+4*t22-2*t24
      t27 = 1/l2
      t29 = Dsquar321(z,l1,l2,0.D0,l4)
      t30 = 1/t29
      t32 = s23r2(z,t12,l2,0.D0,l4)
      t34 = s1vr2(z,t12,l2,0.D0,l4)
      t36 = 1/t32/t34
      t37 = f1(0.D0,z,t12)
      tncoff3 = -t26*t27*t30*t36*t37/1024+t36*t30*t27/l3*t37/512
      t1 = l1*l2
      t2 = f1(l3,z,t1)
      t3 = s1vr2(z,t1,l2,l3,l4)
      t6 = s23r2(z,t1,l2,l3,l4)
      t9 = Dsquar321(z,l1,l2,l3,l4)
      tncoff4 = -t2/t3/t6/t9/l2/l3/512
      t2 = dlog(1-l4)
      t4 = dlog(l1)
      t6 = dlog(l4)
      t8 = dlog(2.D0)
      t10 = f2(0.D0,z)
      t11 = dlog(t10)
      t13 = dlog(z)
      t16 = dlog(1-z)
      t19 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t20 = 1/t19
      t22 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t24 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t27 = f1(0.D0,z,0.D0)
      t28 = 1/t22/t24*t27
      tncoff5 = (t2/2048+t4/2048+t6/2048+t8/512-t11/2048-t13/1024+t16/51
     #2)*t20*t28-1/l3*t20*t28/1024-1/l2*t20*t28/1024
      t2 = Dsquar321(z,l1,0.D0,l3,l4)
      t5 = s23r2(z,0.D0,0.D0,l3,l4)
      t7 = s1vr2(z,0.D0,0.D0,l3,l4)
      t10 = f1(l3,z,0.D0)
      tncoff6 = 1/l3/t2/t5/t7*t10/1024
      t2 = Dsquar321(z,l1,l2,0.D0,l4)
      t5 = l1*l2
      t6 = s23r2(z,t5,l2,0.D0,l4)
      t8 = s1vr2(z,t5,l2,0.D0,l4)
      t11 = f1(0.D0,z,t5)
      tncoff7 = 1/l2/t2/t6/t8*t11/1024
      t1 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t6 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t8 = f1(0.D0,z,0.D0)
      tncoff8 = -1/t1/t3/t6*t8/2048

CCC Numerator functions

      tcnum1 = sq3num3e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num3e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num3e0(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num3e0(var,ndim,z,ux,l1*l2,l2,l3,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num3e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num3e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num3e1(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num3e2(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff8)

CCC Sum

      sq35e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq35e1int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
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
      double precision sq3num3e0
      external sq3num3e0
      double precision sq3num3e1
      external sq3num3e1
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l4 = jac*l4+buf
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

      t2 = dlog(1-l4)
      t4 = dlog(l1)
      t6 = dlog(l4)
      t8 = dlog(2.D0)
      t10 = f2(0.D0,z)
      t11 = dlog(t10)
      t13 = dlog(z)
      t16 = dlog(1-z)
      t19 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t20 = 1/t19
      t22 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t24 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t27 = f1(0.D0,z,0.D0)
      t28 = 1/t22/t24*t27
      tncoff1 = (t2/2048+t4/2048+t6/2048+t8/512-t11/2048-t13/1024+t16/51
     #2)*t20*t28-1/l3*t20*t28/1024-1/l2*t20*t28/1024
      t2 = Dsquar321(z,l1,0.D0,l3,l4)
      t5 = s23r2(z,0.D0,0.D0,l3,l4)
      t7 = s1vr2(z,0.D0,0.D0,l3,l4)
      t10 = f1(l3,z,0.D0)
      tncoff2 = 1/l3/t2/t5/t7*t10/1024
      t2 = Dsquar321(z,l1,l2,0.D0,l4)
      t5 = l1*l2
      t6 = s23r2(z,t5,l2,0.D0,l4)
      t8 = s1vr2(z,t5,l2,0.D0,l4)
      t11 = f1(0.D0,z,t5)
      tncoff3 = 1/l2/t2/t6/t8*t11/1024
      t1 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t6 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t8 = f1(0.D0,z,0.D0)
      tncoff4 = -1/t1/t3/t6*t8/2048

CCC Numerator functions

      tcnum1 = sq3num3e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num3e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num3e0(var,ndim,z,ux,l1*l2,l2,0.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num3e1(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      sq35e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq35e2int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
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
      double precision sq3num3e0
      external sq3num3e0
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l4 = jac*l4+buf
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

      t1 = Dsquar321(z,l1,0.D0,0.D0,l4)
      t3 = s23r2(z,0.D0,0.D0,0.D0,l4)
      t6 = s1vr2(z,0.D0,0.D0,0.D0,l4)
      t8 = f1(0.D0,z,0.D0)
      tncoff1 = -1/t1/t3/t6*t8/2048

CCC Numerator functions

      tcnum1 = sq3num3e0(var,ndim,z,ux,0.D0,0.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      sq35e2int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq35e3int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
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
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l4 = jac*l4+buf
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

      sq35e3int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq35e4int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
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
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l4 = jac*l4+buf
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

      sq35e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

