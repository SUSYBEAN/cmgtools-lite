CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function sq31e0int(var,ndim,Tcoeff)

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
      double precision sq3num1e0
      external sq3num1e0
      double precision sq3num1e1
      external sq3num1e1
      double precision sq3num1e2
      external sq3num1e2
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

      t1 = dlog(l1)
      t3 = dlog(1-z)
      t6 = f3(0.D0,z)
      t8 = dlog(1-t6)
      t9 = f2(0.D0,z)
      t10 = dlog(t9)
      t14 = dlog(1-l4)
      t17 = dlog(2.D0)
      t20 = dlog(z)
      t35 = dlog(l4)
      t50 = t1*t3/512-t8*t10/2048+t14*t1/2048-t10*t17/512-t1*t20/1024-t3
     #*t20/256+t10*t20/1024+t8*t17/512+t8*t3/512-t14*t20/1024-t20*t17/25
     #6+t35*t8/2048+t3*t17/128-t35*t10/2048+t35*t1/2048-t35*t20/1024+t35
     #*t3/512-t14*t10/2048
      t55 = t20**2
      t57 = t10**2
      t71 = t3**2
      t73 = t8**2
      t75 = t14**2
      t77 = t35**2
      t79 = t1**2
      t85 = 0.3141592653589793D1**2
      t87 = t17**2
      t89 = -t1*t10/2048+t14*t8/2048+t55/1024+t57/4096+t1*t17/512-t8*t20
     #/1024+t14*t17/512+t14*t35/2048+t14*t3/512-t3*t10/512+t71/256+t73/4
     #096+t75/4096+t77/4096+t79/4096+t35*t17/512+t1*t8/2048-t85/6144+t87
     #/256
      t91 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t92 = 1/t91
      t94 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t96 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t97 = t96**2
      t99 = 1/t94/t97
      t100 = f1(0.D0,z,1.D0)
      t101 = t99*t100
      t103 = 4*t3
      t104 = 4*t17
      t105 = 2*t20
      t106 = dlog(l3)
      t109 = 1/l3
      t114 = dlog(l2)
      t117 = 1/l2
      tncoff1 = (t50+t89)*t92*t101+(-t1-t103-t8+t10-t104+t105-2*t106-t14
     #-t35)*t109*t92*t101/1024+(-t14-t35-2*t114-t1-t103-t8+t10-t104+t105
     #)*t117*t92*t101/1024+t99*t92*t117*t109*t100/512
      t2 = dlog(1-l3)
      t4 = f3(l3,z)
      t6 = dlog(1-t4)
      t8 = dlog(1-z)
      t10 = dlog(l1)
      t11 = f2(l3,z)
      t12 = dlog(t11)
      t13 = dlog(2.D0)
      t15 = dlog(z)
      t17 = dlog(l3)
      t20 = dlog(1-l4)
      t21 = dlog(l4)
      t23 = 1/l3
      t25 = Dsquar311(z,l1,0.D0,l3,l4)
      t26 = 1/t25
      t28 = s23r2(z,1.D0,1.D0,l3,l4)
      t30 = s1vr2(z,1.D0,1.D0,l3,l4)
      t31 = t30**2
      t33 = 1/t28/t31
      t34 = f1(l3,z,1.D0)
      tncoff2 = (2*t2+t6+4*t8+t10-t12+4*t13-2*t15+2*t17+t20+t21)*t23*t26
     #*t33*t34/1024-t33*t26/l2*t23*t34/512
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = dlog(l2)
      t6 = dlog(l1)
      t7 = f3(0.D0,z)
      t11 = dlog(1-t7+t7*l1*l2)
      t13 = dlog(1-z)
      t15 = 1-l2
      t16 = dlog(t15)
      t17 = f2(0.D0,z)
      t18 = dlog(t17)
      t19 = dlog(2.D0)
      t21 = dlog(z)
      t24 = 1-l1*l2
      t25 = dlog(t24)
      t26 = -t2-t3-2*t4-t6-t11-4*t13-t16+t18-4*t19+2*t21-t25
      t27 = 1/l2
      t29 = Dsquar311(z,l1,l2,0.D0,l4)
      t30 = 1/t29
      t32 = s23r2(z,t24,t15,0.D0,l4)
      t34 = s1vr2(z,t24,t15,0.D0,l4)
      t35 = t34**2
      t37 = 1/t32/t35
      t38 = f1(0.D0,z,t24)
      tncoff3 = -t26*t27*t30*t37*t38/1024-t37*t30*t27/l3*t38/512
      t2 = 1-l1*l2
      t3 = f1(l3,z,t2)
      t4 = 1-l2
      t5 = s1vr2(z,t2,t4,l3,l4)
      t6 = t5**2
      t9 = s23r2(z,t2,t4,l3,l4)
      t12 = Dsquar311(z,l1,l2,l3,l4)
      tncoff4 = t3/t6/t9/t12/l2/l3/512
      t1 = dlog(l1)
      t3 = dlog(2.D0)
      t5 = dlog(z)
      t8 = dlog(1-l4)
      t10 = dlog(l4)
      t12 = f3(0.D0,z)
      t14 = dlog(1-t12)
      t17 = dlog(1-z)
      t19 = f2(0.D0,z)
      t20 = dlog(t19)
      t23 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t24 = 1/t23
      t26 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t28 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t29 = t28**2
      t32 = f1(0.D0,z,1.D0)
      t33 = 1/t26/t29*t32
      tncoff5 = (-t1/2048-t3/512+t5/1024-t8/2048-t10/2048-t14/2048-t17/5
     #12+t20/2048)*t24*t33+1/l3*t24*t33/1024+1/l2*t24*t33/1024
      t2 = Dsquar311(z,l1,0.D0,l3,l4)
      t5 = s23r2(z,1.D0,1.D0,l3,l4)
      t7 = s1vr2(z,1.D0,1.D0,l3,l4)
      t8 = t7**2
      t11 = f1(l3,z,1.D0)
      tncoff6 = -1/l3/t2/t5/t8*t11/1024
      t2 = Dsquar311(z,l1,l2,0.D0,l4)
      t6 = 1-l1*l2
      t7 = 1-l2
      t8 = s23r2(z,t6,t7,0.D0,l4)
      t10 = s1vr2(z,t6,t7,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,t6)
      tncoff7 = -1/l2/t2/t8/t11*t14/1024
      t1 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t3 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t6 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,1.D0)
      tncoff8 = 1/t1/t3/t7*t9/2048

CCC Numerator functions

      tcnum1 = sq3num1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num1e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num1e0(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num1e0(var,ndim,z,ux,1-l1*l2,1-l2,l3,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num1e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num1e1(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num1e1(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num1e2(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff8)

CCC Sum

      sq31e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq31e1int(var,ndim,Tcoeff)

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
      double precision sq3num1e0
      external sq3num1e0
      double precision sq3num1e1
      external sq3num1e1
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

      t1 = dlog(l1)
      t3 = dlog(2.D0)
      t5 = dlog(z)
      t8 = dlog(1-l4)
      t10 = dlog(l4)
      t12 = f3(0.D0,z)
      t14 = dlog(1-t12)
      t17 = dlog(1-z)
      t19 = f2(0.D0,z)
      t20 = dlog(t19)
      t23 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t24 = 1/t23
      t26 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t28 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t29 = t28**2
      t32 = f1(0.D0,z,1.D0)
      t33 = 1/t26/t29*t32
      tncoff1 = (-t1/2048-t3/512+t5/1024-t8/2048-t10/2048-t14/2048-t17/5
     #12+t20/2048)*t24*t33+1/l3*t24*t33/1024+1/l2*t24*t33/1024
      t2 = Dsquar311(z,l1,0.D0,l3,l4)
      t5 = s23r2(z,1.D0,1.D0,l3,l4)
      t7 = s1vr2(z,1.D0,1.D0,l3,l4)
      t8 = t7**2
      t11 = f1(l3,z,1.D0)
      tncoff2 = -1/l3/t2/t5/t8*t11/1024
      t2 = Dsquar311(z,l1,l2,0.D0,l4)
      t6 = 1-l1*l2
      t7 = 1-l2
      t8 = s23r2(z,t6,t7,0.D0,l4)
      t10 = s1vr2(z,t6,t7,0.D0,l4)
      t11 = t10**2
      t14 = f1(0.D0,z,t6)
      tncoff3 = -1/l2/t2/t8/t11*t14/1024
      t1 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t3 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t6 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,1.D0)
      tncoff4 = 1/t1/t3/t7*t9/2048

CCC Numerator functions

      tcnum1 = sq3num1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num1e0(var,ndim,z,ux,1.D0,1.D0,l3,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num1e0(var,ndim,z,ux,1-l1*l2,1-l2,0.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num1e1(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      sq31e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq31e2int(var,ndim,Tcoeff)

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
      double precision sq3num1e0
      external sq3num1e0
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

      t1 = Dsquar311(z,l1,0.D0,0.D0,l4)
      t3 = s23r2(z,1.D0,1.D0,0.D0,l4)
      t6 = s1vr2(z,1.D0,1.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,1.D0)
      tncoff1 = 1/t1/t3/t7*t9/2048

CCC Numerator functions

      tcnum1 = sq3num1e0(var,ndim,z,ux,1.D0,1.D0,0.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      sq31e2int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq31e3int(var,ndim,Tcoeff)

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

      sq31e3int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq31e4int(var,ndim,Tcoeff)

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

      sq31e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

