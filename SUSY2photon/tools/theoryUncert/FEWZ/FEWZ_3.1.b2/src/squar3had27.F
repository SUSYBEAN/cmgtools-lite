CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function sq37e0int(var,ndim,Tcoeff)

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
      double precision sq3num4e0
      external sq3num4e0
      double precision sq3num4e1
      external sq3num4e1
      double precision sq3num4e2
      external sq3num4e2
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
      t3 = t2**2
      t6 = dlog(1-z)
      t7 = t6**2
      t9 = dlog(l4)
      t10 = t9**2
      t12 = dlog(l1)
      t13 = t12**2
      t15 = 0.3141592653589793D1**2
      t17 = dlog(2.D0)
      t18 = t17**2
      t20 = f2(1.D0,z)
      t21 = dlog(t20)
      t28 = t21**2
      t30 = dlog(z)
      t39 = -t3/4096-t7/256-t10/4096-t13/4096+t15/6144-t18/256+t9*t21/20
     #48+t2*t21/2048-t9*t12/2048-t28/4096+t6*t30/256-t6*t17/128-t9*t17/5
     #12-t2*t9/2048
      t64 = t30**2
      t70 = -t2*t17/512-t2*t6/512+t2*t30/1024-t9*t6/512+t12*t30/1024-t12
     #*t17/512+t6*t21/512+t21*t17/512-t21*t30/1024+t30*t17/256+t9*t30/10
     #24-t2*t12/2048-t64/1024-t12*t6/512+t12*t21/2048
      t72 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t73 = 1/t72
      t75 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t77 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t79 = 1/t75/t77
      t80 = f1(1.D0,z,0.D0)
      t81 = t79*t80
      t83 = 1-l3
      t84 = dlog(t83)
      t86 = 4*t6
      t87 = 4*t17
      t88 = 2*t30
      t91 = -1/t83
      t96 = dlog(l2)
      t99 = 1/l2
      tncoff1 = (t39+t70)*t73*t81+(-t9-t2-2*t84-t12-t86+t21-t87+t88)*t91
     #*t73*t81/1024-(-t2-t9-2*t96-t12-t86+t21-t87+t88)*t99*t73*t81/1024+
     #t79*t73*t99*t91*t80/512
      t1 = dlog(l4)
      t3 = dlog(1-l4)
      t4 = 1-l3
      t5 = dlog(t4)
      t7 = dlog(l1)
      t9 = dlog(1-z)
      t11 = dlog(l3)
      t13 = dlog(2.D0)
      t15 = dlog(z)
      t17 = f2(l3,z)
      t18 = dlog(t17)
      t21 = -1/t4
      t23 = Dsquar321(z,l1,0.D0,l3,l4)
      t24 = 1/t23
      t26 = s2vr2(z,0.D0,0.D0,l3,l4)
      t28 = s23r2(z,0.D0,0.D0,l3,l4)
      t30 = 1/t26/t28
      t31 = f1(l3,z,0.D0)
      tncoff2 = -(-t1-t3-2*t5-t7-4*t9-2*t11-4*t13+2*t15+t18)*t21*t24*t30
     #*t31/1024-t30*t24/l2*t21*t31/512
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = dlog(l2)
      t6 = dlog(l1)
      t8 = dlog(1-z)
      t11 = dlog(1-l2)
      t12 = l1*l2
      t14 = dlog(1-t12)
      t15 = f3(1.D0,z)
      t19 = dlog(1-t15*l1*l2)
      t20 = dlog(2.D0)
      t22 = dlog(z)
      t24 = f2(1.D0,z)
      t25 = dlog(t24)
      t26 = -t2-t3-2*t4-t6-4*t8-t11-t14-t19-4*t20+2*t22+t25
      t27 = 1/l2
      t29 = Dsquar321(z,l1,l2,1.D0,l4)
      t30 = 1/t29
      t32 = s2vr2(z,t12,l2,1.D0,l4)
      t34 = s23r2(z,t12,l2,1.D0,l4)
      t36 = 1/t32/t34
      t37 = f1(1.D0,z,t12)
      tncoff3 = t26*t27*t30*t36*t37/1024-t36*t30*t27/(l3-1)*t37/512
      t1 = l1*l2
      t2 = f1(l3,z,t1)
      t3 = s2vr2(z,t1,l2,l3,l4)
      t6 = s23r2(z,t1,l2,l3,l4)
      t9 = Dsquar321(z,l1,l2,l3,l4)
      tncoff4 = t2/t3/t6/t9/l2/(l3-1)/512
      t2 = dlog(1-l4)
      t4 = dlog(l4)
      t7 = dlog(1-z)
      t9 = f2(1.D0,z)
      t10 = dlog(t9)
      t12 = dlog(z)
      t14 = dlog(2.D0)
      t16 = dlog(l1)
      t19 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t20 = 1/t19
      t22 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t24 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t27 = f1(1.D0,z,0.D0)
      t28 = 1/t22/t24*t27
      tncoff5 = (t2/2048+t4/2048+t7/512-t10/2048-t12/1024+t14/512+t16/20
     #48)*t20*t28+1/(l3-1)*t20*t28/1024-1/l2*t20*t28/1024
      t3 = Dsquar321(z,l1,0.D0,l3,l4)
      t6 = s2vr2(z,0.D0,0.D0,l3,l4)
      t8 = s23r2(z,0.D0,0.D0,l3,l4)
      t11 = f1(l3,z,0.D0)
      tncoff6 = -1/(l3-1)/t3/t6/t8*t11/1024
      t2 = Dsquar321(z,l1,l2,1.D0,l4)
      t5 = l1*l2
      t6 = s2vr2(z,t5,l2,1.D0,l4)
      t8 = s23r2(z,t5,l2,1.D0,l4)
      t11 = f1(1.D0,z,t5)
      tncoff7 = 1/l2/t2/t6/t8*t11/1024
      t1 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t3 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t6 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t8 = f1(1.D0,z,0.D0)
      tncoff8 = -1/t1/t3/t6*t8/2048

CCC Numerator functions

      tcnum1 = sq3num4e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num4e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num4e0(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num4e0(var,ndim,z,ux,l1*l2,l2,l3,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num4e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num4e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num4e1(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num4e2(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff8)

CCC Sum

      sq37e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq37e1int(var,ndim,Tcoeff)

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
      double precision sq3num4e0
      external sq3num4e0
      double precision sq3num4e1
      external sq3num4e1
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
      t4 = dlog(l4)
      t7 = dlog(1-z)
      t9 = f2(1.D0,z)
      t10 = dlog(t9)
      t12 = dlog(z)
      t14 = dlog(2.D0)
      t16 = dlog(l1)
      t19 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t20 = 1/t19
      t22 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t24 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t27 = f1(1.D0,z,0.D0)
      t28 = 1/t22/t24*t27
      tncoff1 = (t2/2048+t4/2048+t7/512-t10/2048-t12/1024+t14/512+t16/20
     #48)*t20*t28+1/(l3-1)*t20*t28/1024-1/l2*t20*t28/1024
      t3 = Dsquar321(z,l1,0.D0,l3,l4)
      t6 = s2vr2(z,0.D0,0.D0,l3,l4)
      t8 = s23r2(z,0.D0,0.D0,l3,l4)
      t11 = f1(l3,z,0.D0)
      tncoff2 = -1/(l3-1)/t3/t6/t8*t11/1024
      t2 = Dsquar321(z,l1,l2,1.D0,l4)
      t5 = l1*l2
      t6 = s2vr2(z,t5,l2,1.D0,l4)
      t8 = s23r2(z,t5,l2,1.D0,l4)
      t11 = f1(1.D0,z,t5)
      tncoff3 = 1/l2/t2/t6/t8*t11/1024
      t1 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t3 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t6 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t8 = f1(1.D0,z,0.D0)
      tncoff4 = -1/t1/t3/t6*t8/2048

CCC Numerator functions

      tcnum1 = sq3num4e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num4e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num4e0(var,ndim,z,ux,l1*l2,l2,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num4e1(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      sq37e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq37e2int(var,ndim,Tcoeff)

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
      double precision sq3num4e0
      external sq3num4e0
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

      t1 = Dsquar321(z,l1,0.D0,1.D0,l4)
      t3 = s2vr2(z,0.D0,0.D0,1.D0,l4)
      t6 = s23r2(z,0.D0,0.D0,1.D0,l4)
      t8 = f1(1.D0,z,0.D0)
      tncoff1 = -1/t1/t3/t6*t8/2048

CCC Numerator functions

      tcnum1 = sq3num4e0(var,ndim,z,ux,0.D0,0.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      sq37e2int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq37e3int(var,ndim,Tcoeff)

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

      sq37e3int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq37e4int(var,ndim,Tcoeff)

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

      sq37e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

