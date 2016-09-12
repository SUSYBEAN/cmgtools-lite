CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function sq33e0int(var,ndim,Tcoeff)

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
      double precision sq3num2e0
      external sq3num2e0
      double precision sq3num2e1
      external sq3num2e1
      double precision sq3num2e2
      external sq3num2e2
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
      t6 = f3(1.D0,z)
      t8 = dlog(1-t6)
      t9 = t8**2
      t11 = f2(1.D0,z)
      t12 = dlog(t11)
      t13 = t12**2
      t15 = dlog(z)
      t16 = t15**2
      t19 = dlog(1-z)
      t20 = t19**2
      t23 = dlog(1-l4)
      t24 = t23**2
      t26 = dlog(l4)
      t27 = t26**2
      t29 = dlog(l1)
      t30 = t29**2
      t50 = -t1/6144+t4/256+t9/4096+t13/4096+t16/1024+t20/256+t24/4096+t
     #27/4096+t30/4096+t23*t8/2048+t23*t29/2048-t15*t3/256-t23*t15/1024-
     #t19*t15/256+t26*t3/512+t26*t19/512-t26*t15/1024+t19*t3/128
      t89 = -t19*t12/512-t12*t3/512+t23*t19/512+t23*t26/2048+t23*t3/512+
     #t8*t19/512-t23*t12/2048+t29*t19/512+t8*t3/512-t8*t12/2048+t12*t15/
     #1024-t8*t15/1024-t29*t15/1024+t29*t3/512-t29*t12/2048+t29*t8/2048-
     #t26*t12/2048+t26*t8/2048+t26*t29/2048
      t91 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t92 = 1/t91
      t94 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t95 = t94**2
      t97 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t99 = 1/t95/t97
      t100 = f1(1.D0,z,1.D0)
      t101 = t99*t100
      t103 = 4*t19
      t104 = 4*t3
      t105 = 2*t15
      t106 = dlog(l3)
      t109 = 1/l3
      t114 = dlog(l2)
      t117 = 1/l2
      tncoff1 = (t50+t89)*t92*t101-(t26+t23+t29+t8+t103-t12+t104-t105+2*
     #t106)*t109*t92*t101/1024-(t23+t26+2*t114+t29+t8+t103-t12+t104-t105
     #)*t117*t92*t101/1024+t99*t92*t117*t109*t100/512
      t1 = dlog(l4)
      t3 = dlog(1-l4)
      t4 = 1-l3
      t5 = f3(t4,z)
      t7 = dlog(1-t5)
      t9 = dlog(1-z)
      t11 = dlog(l1)
      t12 = f2(t4,z)
      t13 = dlog(t12)
      t14 = dlog(2.D0)
      t16 = dlog(z)
      t18 = dlog(t4)
      t20 = dlog(l3)
      t23 = 1/l3
      t25 = Dsquar311(z,l1,0.D0,t4,l4)
      t26 = 1/t25
      t28 = s2vr2(z,1.D0,1.D0,t4,l4)
      t29 = t28**2
      t31 = s23r2(z,1.D0,1.D0,t4,l4)
      t33 = 1/t29/t31
      t34 = f1(t4,z,1.D0)
      tncoff2 = (t1+t3+t7+4*t9+t11-t13+4*t14-2*t16+2*t18+2*t20)*t23*t26*
     #t33*t34/1024-t33*t26/l2*t23*t34/512
      t2 = dlog(1-l4)
      t3 = dlog(l4)
      t4 = dlog(l2)
      t6 = dlog(l1)
      t7 = f3(1.D0,z)
      t11 = dlog(1-t7+t7*l1*l2)
      t13 = dlog(1-z)
      t15 = 1-l2
      t16 = dlog(t15)
      t17 = f2(1.D0,z)
      t18 = dlog(t17)
      t19 = dlog(2.D0)
      t21 = dlog(z)
      t24 = 1-l1*l2
      t25 = dlog(t24)
      t26 = t2+t3+2*t4+t6+t11+4*t13+t16-t18+4*t19-2*t21+t25
      t27 = 1/l2
      t29 = Dsquar311(z,l1,l2,1.D0,l4)
      t30 = 1/t29
      t32 = s2vr2(z,t24,t15,1.D0,l4)
      t33 = t32**2
      t34 = 1/t33
      t35 = s23r2(z,t24,t15,1.D0,l4)
      t36 = 1/t35
      t38 = f1(1.D0,z,t24)
      tncoff3 = t26*t27*t30*t34*t36*t38/1024-t27*t34*t36*t30/l3*t38/512
      t1 = 1-l3
      t3 = 1-l1*l2
      t4 = f1(t1,z,t3)
      t5 = 1-l2
      t6 = s2vr2(z,t3,t5,t1,l4)
      t7 = t6**2
      t10 = s23r2(z,t3,t5,t1,l4)
      t13 = Dsquar311(z,l1,l2,t1,l4)
      tncoff4 = t4/t7/t10/t13/l2/l3/512
      t2 = dlog(1-z)
      t4 = f3(1.D0,z)
      t6 = dlog(1-t4)
      t8 = dlog(l1)
      t11 = dlog(1-l4)
      t13 = dlog(l4)
      t15 = dlog(2.D0)
      t17 = dlog(z)
      t19 = f2(1.D0,z)
      t20 = dlog(t19)
      t23 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t24 = 1/t23
      t26 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t27 = t26**2
      t29 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t32 = f1(1.D0,z,1.D0)
      t33 = 1/t27/t29*t32
      tncoff5 = (-t2/512-t6/2048-t8/2048-t11/2048-t13/2048-t15/512+t17/1
     #024+t20/2048)*t24*t33+1/l3*t24*t33/1024+1/l2*t24*t33/1024
      t2 = 1-l3
      t3 = Dsquar311(z,l1,0.D0,t2,l4)
      t6 = s2vr2(z,1.D0,1.D0,t2,l4)
      t7 = t6**2
      t9 = s23r2(z,1.D0,1.D0,t2,l4)
      t12 = f1(t2,z,1.D0)
      tncoff6 = -1/l3/t3/t7/t9*t12/1024
      t2 = Dsquar311(z,l1,l2,1.D0,l4)
      t6 = 1-l1*l2
      t7 = 1-l2
      t8 = s2vr2(z,t6,t7,1.D0,l4)
      t9 = t8**2
      t11 = s23r2(z,t6,t7,1.D0,l4)
      t14 = f1(1.D0,z,t6)
      tncoff7 = -1/l2/t2/t9/t11*t14/1024
      t1 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t3 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t4 = t3**2
      t7 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t9 = f1(1.D0,z,1.D0)
      tncoff8 = 1/t1/t4/t7*t9/2048

CCC Numerator functions

      tcnum1 = sq3num2e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num2e0(var,ndim,z,ux,1.D0,1.D0,1-l3,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num2e0(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num2e0(var,ndim,z,ux,1-l1*l2,1-l2,1-l3,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = sq3num2e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = sq3num2e1(var,ndim,z,ux,1.D0,1.D0,1-l3,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = sq3num2e1(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = sq3num2e2(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff8)

CCC Sum

      sq33e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq33e1int(var,ndim,Tcoeff)

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
      double precision sq3num2e0
      external sq3num2e0
      double precision sq3num2e1
      external sq3num2e1
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

      t2 = dlog(1-z)
      t4 = f3(1.D0,z)
      t6 = dlog(1-t4)
      t8 = dlog(l1)
      t11 = dlog(1-l4)
      t13 = dlog(l4)
      t15 = dlog(2.D0)
      t17 = dlog(z)
      t19 = f2(1.D0,z)
      t20 = dlog(t19)
      t23 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t24 = 1/t23
      t26 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t27 = t26**2
      t29 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t32 = f1(1.D0,z,1.D0)
      t33 = 1/t27/t29*t32
      tncoff1 = (-t2/512-t6/2048-t8/2048-t11/2048-t13/2048-t15/512+t17/1
     #024+t20/2048)*t24*t33+1/l3*t24*t33/1024+1/l2*t24*t33/1024
      t2 = 1-l3
      t3 = Dsquar311(z,l1,0.D0,t2,l4)
      t6 = s2vr2(z,1.D0,1.D0,t2,l4)
      t7 = t6**2
      t9 = s23r2(z,1.D0,1.D0,t2,l4)
      t12 = f1(t2,z,1.D0)
      tncoff2 = -1/l3/t3/t7/t9*t12/1024
      t2 = Dsquar311(z,l1,l2,1.D0,l4)
      t6 = 1-l1*l2
      t7 = 1-l2
      t8 = s2vr2(z,t6,t7,1.D0,l4)
      t9 = t8**2
      t11 = s23r2(z,t6,t7,1.D0,l4)
      t14 = f1(1.D0,z,t6)
      tncoff3 = -1/l2/t2/t9/t11*t14/1024
      t1 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t3 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t4 = t3**2
      t7 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t9 = f1(1.D0,z,1.D0)
      tncoff4 = 1/t1/t4/t7*t9/2048

CCC Numerator functions

      tcnum1 = sq3num2e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = sq3num2e0(var,ndim,z,ux,1.D0,1.D0,1-l3,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = sq3num2e0(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,l8,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = sq3num2e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      sq33e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq33e2int(var,ndim,Tcoeff)

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
      double precision sq3num2e0
      external sq3num2e0
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

      t1 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t3 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t4 = t3**2
      t7 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t9 = f1(1.D0,z,1.D0)
      tncoff1 = 1/t1/t4/t7*t9/2048

CCC Numerator functions

      tcnum1 = sq3num2e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      sq33e2int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq33e3int(var,ndim,Tcoeff)

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

      sq33e3int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function sq33e4int(var,ndim,Tcoeff)

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

      sq33e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

