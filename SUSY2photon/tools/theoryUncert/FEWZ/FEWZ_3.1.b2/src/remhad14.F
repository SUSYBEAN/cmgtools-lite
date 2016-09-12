CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function rm14e0int(var,ndim,Tcoeff)

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
      double precision rm1num3e0
      external rm1num3e0
      double precision rm1num3e1
      external rm1num3e1
      double precision rm1num3e2
      external rm1num3e2
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

      t1 = -1+z
      t2 = dlog(l4)
      t3 = t2**2
      t6 = dlog(1-l4)
      t7 = t6**2
      t9 = dlog(l1)
      t10 = t9**2
      t13 = dlog(1-l1)
      t14 = t13**2
      t22 = dlog(2.D0)
      t31 = -3*t3-3*t7-3*t10-3*t14-6*t9*t13-6*t2*t13-6*t6*t13-18*t13*t22
     #-6*t6*t2-18*t2*t22-18*t6*t22
      t39 = dlog(-t1)
      t46 = f2(0.D0,z)
      t47 = dlog(t46)
      t52 = dlog(z)
      t59 = f3(0.D0,z)
      t62 = dlog(1-l1*t59)
      t65 = -18*t9*t22-6*t2*t9-6*t6*t9-72*t39*t22-24*t6*t39-24*t2*t39+6*
     #t9*t47-24*t9*t39+12*t9*t52+12*t13*t52-24*t39*t13-18*t62*t22
      t89 = 12*t62*t52+6*t62*t47-24*t62*t39-6*t62*t13+6*t47*t13-6*t9*t62
     #-6*t2*t62-6*t6*t62+6*t6*t47+12*t6*t52+6*t2*t47
      t102 = t22**2
      t104 = t39**2
      t106 = t47**2
      t108 = t52**2
      t110 = t62**2
      t112 = 0.3141592653589793D1**2
      t114 = 12*t2*t52+24*t39*t47+48*t39*t52+36*t52*t22+18*t47*t22-12*t4
     #7*t52-27*t102-48*t104-3*t106-12*t108-3*t110+2*t112
      t118 = s23r2(z,l1,0.D0,0.D0,l4)
      t119 = t118**2
      t121 = s1vr2(z,l1,0.D0,0.D0,l4)
      t122 = t121**2
      t125 = f1(0.D0,z,l1)
      t126 = 1/t119/t122*t125
      t129 = dlog(l3)
      t131 = 4*t39
      t132 = 2*t52
      t133 = 3*t22
      t136 = 1/l3
      t140 = dlog(l2)
      t143 = 1/l2
      tncoff1 = t1*(t31+t65+t89+t114)*t126/6144-t1*(-2*t129-t62+t47-t131
     #-t13-t9+t132-t133-t6-t2)*t136*t126/512-(-t6-t2-t140+t47-t62-t131-t
     #133-t13-t9+t132)*t1*t143*t126/1024-t1*t143*t136*t126/512
      t1 = -1+z
      t2 = l3-1
      t3 = t1*t2
      t4 = dlog(l3)
      t6 = f2(l3,z)
      t7 = dlog(t6)
      t9 = dlog(-t2)
      t11 = f3(l3,z)
      t14 = dlog(1-l1*t11)
      t16 = dlog(-t1)
      t18 = dlog(z)
      t21 = dlog(1-l1)
      t23 = dlog(1-l4)
      t24 = dlog(l1)
      t25 = dlog(l4)
      t26 = dlog(2.D0)
      t28 = 2*t4-t7+2*t9+t14+4*t16-2*t18+t21+t23+t24+t25+3*t26
      t31 = s23r2(z,l1,0.D0,l3,l4)
      t32 = t31**2
      t35 = s1vr2(z,l1,0.D0,l3,l4)
      t36 = t35**2
      t38 = f1(l3,z,l1)
      t40 = 1/l3/t32/t36*t38
      tncoff2 = t3*t28*t40/512-t3/l2*t40/512
      t1 = -1+z
      t3 = dlog(1-l4)
      t4 = dlog(l4)
      t5 = dlog(l2)
      t7 = dlog(1-l1)
      t8 = dlog(l1)
      t9 = 2-l2
      t10 = dlog(t9)
      t11 = dlog(z)
      t13 = f3(0.D0,z)
      t16 = dlog(1-l1*t13)
      t17 = dlog(2.D0)
      t20 = dlog(-t1)
      t22 = f2(0.D0,z)
      t23 = dlog(t22)
      t24 = t3+t4+t5+t7+t8+t10-2*t11+t16+2*t17+4*t20-t23
      t27 = -1/t9
      t29 = 1/l2
      t30 = l2/2
      t31 = s23r2(z,l1,t30,0.D0,l4)
      t32 = t31**2
      t33 = 1/t32
      t35 = s1vr2(z,l1,t30,0.D0,l4)
      t36 = t35**2
      t37 = 1/t36
      t38 = f1(0.D0,z,l1)
      tncoff3 = t1*t24*t27*t29*t33*t37*t38/512-t27*t37*t33*t29*t1/l3*t38
     #/256
      t10 = l2/2
      t11 = s23r2(z,l1,t10,l3,l4)
      t12 = t11**2
      t15 = s1vr2(z,l1,t10,l3,l4)
      t16 = t15**2
      t18 = f1(l3,z,l1)
      tncoff4 = -(-1+z)*(l3-1)/(-2+l2)/l2/l3/t12/t16*t18/256
      t1 = -1+z
      t2 = f2(0.D0,z)
      t3 = dlog(t2)
      t6 = dlog(-t1)
      t8 = dlog(2.D0)
      t11 = dlog(1-l1)
      t13 = dlog(l1)
      t15 = dlog(z)
      t17 = f3(0.D0,z)
      t20 = dlog(1-l1*t17)
      t22 = dlog(l4)
      t25 = dlog(1-l4)
      t29 = s23r2(z,l1,0.D0,0.D0,l4)
      t30 = t29**2
      t32 = s1vr2(z,l1,0.D0,0.D0,l4)
      t33 = t32**2
      t36 = f1(0.D0,z,l1)
      t37 = 1/t30/t33*t36
      tncoff5 = t1*(-6*t3+24*t6+18*t8+6*t11+6*t13-12*t15+6*t20+6*t22+6*t
     #25)*t37/6144-t1/l3*t37/512-t1/l2*t37/1024
      t6 = s23r2(z,l1,0.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,0.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff6 = -(-1+z)*(l3-1)/l3/t7/t10*t13/512
      t7 = l2/2
      t8 = s23r2(z,l1,t7,0.D0,l4)
      t9 = t8**2
      t11 = s1vr2(z,l1,t7,0.D0,l4)
      t12 = t11**2
      t15 = f1(0.D0,z,l1)
      tncoff7 = -(-1+z)/(-2+l2)/l2/t9/t12*t15/512
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff8 = -(-1+z)/t3/t7*t9/1024

CCC Numerator functions

      tcnum1 = rm1num3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = rm1num3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = rm1num3e0(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = rm1num3e0(var,ndim,z,ux,l1,l2/2,l3,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = rm1num3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = rm1num3e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = rm1num3e1(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = rm1num3e2(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff8)

CCC Sum

      rm14e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm14e1int(var,ndim,Tcoeff)

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
      double precision rm1num3e0
      external rm1num3e0
      double precision rm1num3e1
      external rm1num3e1
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

      t1 = -1+z
      t2 = dlog(l4)
      t3 = f2(0.D0,z)
      t4 = dlog(t3)
      t6 = dlog(-t1)
      t9 = dlog(1-l1)
      t10 = f3(0.D0,z)
      t13 = dlog(1-l1*t10)
      t14 = dlog(z)
      t16 = dlog(2.D0)
      t18 = dlog(l1)
      t20 = dlog(1-l4)
      t23 = s23r2(z,l1,0.D0,0.D0,l4)
      t24 = t23**2
      t26 = s1vr2(z,l1,0.D0,0.D0,l4)
      t27 = t26**2
      t30 = f1(0.D0,z,l1)
      t31 = 1/t24/t27*t30
      tncoff1 = -t1*(-t2+t4-4*t6-t9-t13+2*t14-3*t16-t18-t20)*t31/1024-t1
     #/l3*t31/512-t1/l2*t31/1024
      t6 = s23r2(z,l1,0.D0,l3,l4)
      t7 = t6**2
      t9 = s1vr2(z,l1,0.D0,l3,l4)
      t10 = t9**2
      t13 = f1(l3,z,l1)
      tncoff2 = -(-1+z)*(l3-1)/l3/t7/t10*t13/512
      t7 = l2/2
      t8 = s23r2(z,l1,t7,0.D0,l4)
      t9 = t8**2
      t11 = s1vr2(z,l1,t7,0.D0,l4)
      t12 = t11**2
      t15 = f1(0.D0,z,l1)
      tncoff3 = -(-1+z)/(-2+l2)/l2/t9/t12*t15/512
      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff4 = -(-1+z)/t3/t7*t9/1024

CCC Numerator functions

      tcnum1 = rm1num3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = rm1num3e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999.D0
     #,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = rm1num3e0(var,ndim,z,ux,l1,l2/2,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = rm1num3e1(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      rm14e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm14e2int(var,ndim,Tcoeff)

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
      double precision rm1num3e0
      external rm1num3e0
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

      t2 = s23r2(z,l1,0.D0,0.D0,l4)
      t3 = t2**2
      t6 = s1vr2(z,l1,0.D0,0.D0,l4)
      t7 = t6**2
      t9 = f1(0.D0,z,l1)
      tncoff1 = -(-1+z)/t3/t7*t9/1024

CCC Numerator functions

      tcnum1 = rm1num3e0(var,ndim,z,ux,l1,0.D0,0.D0,l4,l5,l6,l7,l8,-999.
     #D0,-999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      rm14e2int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm14e3int(var,ndim,Tcoeff)

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

      rm14e3int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm14e4int(var,ndim,Tcoeff)

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

      rm14e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

