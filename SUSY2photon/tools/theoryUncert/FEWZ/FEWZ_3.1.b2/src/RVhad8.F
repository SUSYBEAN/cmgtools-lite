CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function RV8e0int(var,ndim,Tcoeff)

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
      double precision RVnum5e0
      external RVnum5e0
      double precision RVnum5e1
      external RVnum5e1
      double precision RVnum5e2
      external RVnum5e2
      double precision RVnum5e3
      external RVnum5e3
      double precision RVnum5e4
      external RVnum5e4
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

      t2 = 0.3141592653589793D1**2
      t3 = t2**2
      t6 = 1-z
      t7 = dlog(t6)
      t11 = 4*t2
      t12 = t7**2
      t15 = 24*t7
      t16 = 12*t12
      t19 = -1/t6
      t22 = dlog(l1)
      t26 = t22**2
      t32 = 1/l1
      tncoff1 = -z3/48+1.D0/16.D0+t3/2560-t2/96+(-4*t7*t2+8*z3+t11+4*t12
     #*t7+t15-t16-24)*t19/384-(-4*t22*t2+24*t22+t26*t22-6*t26+16*z3+8*t2
     #-48)*t32/1536-(-12*t22+3*t26+12*t7*t22-t11-t15+t16+24)*t32*t19/768
      t1 = dlog(l1)
      t2 = 0.3141592653589793D1**2
      t6 = dlog(1-l1)
      t7 = t6**2
      t10 = t1**2
      t22 = t1*t6
      t25 = 4*t1*t2-3*t1*t7-3*t10*t6+4*t6*t2-16*z3-t7*t6-t10*t1-8*t2-24*
     #t6-24*t1+6*t10+12*t22+6*t7+48
      t26 = 1/l1
      t34 = 1-z
      t35 = dlog(t34)
      t42 = t35**2
      t44 = 4*t2+12*t6+12*t1-3*t10-6*t22-12*t35*t6-12*t35*t1-3*t7+24*t35
     #-12*t42-24
      tncoff2 = -t25*t26/1536+t44*t26/t34/768
      t1 = 1-z
      t2 = dlog(t1)
      t3 = 0.3141592653589793D1**2
      t7 = dlog(z)
      t8 = t7**2
      t11 = t2**2
      t14 = 24*t2
      t15 = 12*t8
      t16 = 4*t3
      t17 = 24*t7
      t18 = 12*t11
      t20 = 24*t7*t2
      t27 = 4*t2*t3-8*z3+4*t8*t7-4*t11*t2-t14+t15-t16+t17+t18-t20-4*t3*t
     #7+12*t11*t7-12*t2*t8+24
      t29 = -1/t1
      t32 = dlog(l1)
      t33 = t32**2
      t40 = t14-t15+t16-t17-t18+t20-3*t33+12*t32*t7-12*t2*t32+12*t32-24
      tncoff3 = t27*t29/384-t40/l1*t29/768
      t1 = dlog(l1)
      t2 = dlog(z)
      t6 = t1**2
      t8 = t2**2
      t11 = dlog(1-l1)
      t12 = t11**2
      t14 = 0.3141592653589793D1**2
      t17 = 1-z
      t18 = dlog(t17)
      t20 = t18**2
      t33 = -12*t1*t2+24-12*t1+3*t6+12*t8+3*t12-4*t14+24*t2-24*t18+12*t2
     #0-12*t11+6*t1*t11-12*t11*t2+12*t18*t11-24*t18*t2+12*t18*t1
      tncoff4 = t33/l1/t17/768
      t1 = 0.3141592653589793D1**2
      t5 = 1-z
      t6 = dlog(t5)
      t7 = t6**2
      t9 = 12*t6
      t12 = -1/t5
      t16 = dlog(l1)
      t17 = t16**2
      t21 = 1/l1
      tncoff5 = 1.D0/32.D0-t1/192-z3/96+(2*t1-6*t7+t9-12)*t12/384-(4*t1-
     #3*t17+12*t16-24)*t21/1536-(-t9-6*t16+12)*t21*t12/768
      t2 = dlog(1-l1)
      t3 = t2**2
      t5 = dlog(l1)
      t8 = 0.3141592653589793D1**2
      t10 = t5**2
      t15 = 1/l1
      t20 = 1-z
      t21 = dlog(t20)
      tncoff6 = -(3*t3+6*t5*t2-4*t8+3*t10-12*t2-12*t5+24)*t15/1536+(6*t2
     #+6*t5+12*t21-12)*t15/t20/768
      t1 = 1-z
      t2 = dlog(t1)
      t3 = dlog(z)
      t6 = t3**2
      t8 = t2**2
      t10 = 0.3141592653589793D1**2
      t12 = 12*t2
      t13 = 12*t3
      t16 = -1/t1
      t19 = dlog(l1)
      tncoff7 = (-12*t2*t3+6*t6+6*t8-2*t10-t12+t13+12)*t16/384-(t12-t13+
     #6*t19-12)/l1*t16/768
      t1 = dlog(l1)
      t3 = dlog(z)
      t6 = dlog(1-l1)
      t8 = 1-z
      t9 = dlog(t8)
      tncoff8 = (-6*t1+12*t3-6*t6-12*t9+12)/l1/t8/768
      t1 = 0.3141592653589793D1**2
      t3 = 1-z
      t4 = dlog(t3)
      t8 = -1/t3
      t11 = dlog(l1)
      t14 = 1/l1
      tncoff9 = -t1/384+1.D0/64.D0+(-6+6*t4)*t8/384-(-12+6*t11)*t14/1536
     #-t14*t8/128
      t1 = dlog(l1)
      t4 = dlog(1-l1)
      t7 = 1/l1
      tncoff10 = -(-6*t1+12-6*t4)*t7/1536+t7/(-1+z)/128
      t1 = dlog(z)
      t2 = 1-z
      t3 = dlog(t2)
      t7 = -1/t2
      tncoff11 = (6*t1-6*t3+6)*t7/384+1/l1*t7/128
      tncoff12 = -1/l1/(-1+z)/128
      tncoff13 = 1.D0/128.D0-1/(-1+z)/128+1/l1/256
      tncoff14 = -1/l1/256
      tncoff15 = 1/(-1+z)/128
      tncoff16 = 1.D0/256.D0

CCC Numerator functions

      tcnum1 = RVnum5e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum5e0(var,ndim,1.D0,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum5e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum5e0(var,ndim,z,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum5e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum5e1(var,ndim,1.D0,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum5e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum5e1(var,ndim,z,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = RVnum5e2(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = RVnum5e2(var,ndim,1.D0,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = RVnum5e2(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = RVnum5e2(var,ndim,z,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = RVnum5e3(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = RVnum5e3(var,ndim,1.D0,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = RVnum5e3(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = RVnum5e4(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff16)

CCC Sum

      RV8e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+t
     #cnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff16

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV8e1int(var,ndim,Tcoeff)

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
      double precision RVnum5e0
      external RVnum5e0
      double precision RVnum5e1
      external RVnum5e1
      double precision RVnum5e2
      external RVnum5e2
      double precision RVnum5e3
      external RVnum5e3
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
      t4 = 1-z
      t5 = dlog(t4)
      t8 = t5**2
      t12 = -1/t4
      t16 = dlog(l1)
      t18 = t16**2
      t21 = 1/l1
      tncoff1 = 1.D0/32.D0-t1/192-z3/96-(-12*t5+12-2*t1+6*t8)*t12/384+(-
     #4*t1+24-12*t16+3*t18)*t21/1536+(2*t5-2+t16)*t21*t12/128
      t2 = dlog(1-l1)
      t4 = 0.3141592653589793D1**2
      t6 = dlog(l1)
      t8 = t6**2
      t10 = t2**2
      t15 = 1/l1
      t18 = 1-z
      t19 = dlog(t18)
      tncoff2 = (-24+12*t2+4*t4+12*t6-3*t8-3*t10-6*t6*t2)*t15/1536-(-2*t
     #19-t2+2-t6)*t15/t18/128
      t1 = 1-z
      t2 = dlog(t1)
      t4 = dlog(z)
      t7 = t4**2
      t9 = 0.3141592653589793D1**2
      t12 = t2**2
      t16 = -1/t1
      t20 = dlog(l1)
      tncoff3 = -(12*t2+12*t2*t4-12-6*t7+2*t9-12*t4-6*t12)*t16/384+(-2*t
     #2+2-t20+2*t4)/l1*t16/128
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
      tncoff5 = -t1/384+1.D0/64.D0-(-6*t4+6)*t8/384+(-6*t11+12)*t14/1536
     #-t14*t8/128
      t2 = dlog(1-l1)
      t4 = dlog(l1)
      t7 = 1/l1
      tncoff6 = (-12+6*t2+6*t4)*t7/1536+t7/(-1+z)/128
      t1 = 1-z
      t2 = dlog(t1)
      t3 = dlog(z)
      t7 = -1/t1
      tncoff7 = -(6*t2-6*t3-6)*t7/384+1/l1*t7/128
      tncoff8 = -1/l1/(-1+z)/128
      tncoff9 = 1.D0/128.D0-1/(-1+z)/128+1/l1/256
      tncoff10 = -1/l1/256
      tncoff11 = 1/(-1+z)/128
      tncoff12 = 1.D0/256.D0

CCC Numerator functions

      tcnum1 = RVnum5e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum5e0(var,ndim,1.D0,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum5e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum5e0(var,ndim,z,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum5e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum5e1(var,ndim,1.D0,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum5e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum5e1(var,ndim,z,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = RVnum5e2(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = RVnum5e2(var,ndim,1.D0,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = RVnum5e2(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = RVnum5e3(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      RV8e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV8e2int(var,ndim,Tcoeff)

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
      double precision RVnum5e0
      external RVnum5e0
      double precision RVnum5e1
      external RVnum5e1
      double precision RVnum5e2
      external RVnum5e2
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
      tncoff1 = 1.D0/64.D0-t1/384+(-2+2*t4)*t8/128-(-2+t11)*t13/256-t13*
     #t8/128
      t2 = dlog(1-l1)
      t3 = dlog(l1)
      t5 = 1/l1
      tncoff2 = -(-t2+2-t3)*t5/256+t5/(-1+z)/128
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

      tcnum1 = RVnum5e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum5e0(var,ndim,1.D0,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum5e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum5e0(var,ndim,z,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = RVnum5e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = RVnum5e1(var,ndim,1.D0,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = RVnum5e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = RVnum5e2(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff8)

CCC Sum

      RV8e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV8e3int(var,ndim,Tcoeff)

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
      double precision RVnum5e0
      external RVnum5e0
      double precision RVnum5e1
      external RVnum5e1
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

      tcnum1 = RVnum5e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = RVnum5e0(var,ndim,1.D0,ux,1-l1,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = RVnum5e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0,9
     #99.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = RVnum5e1(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      RV8e3int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function RV8e4int(var,ndim,Tcoeff)

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
      double precision RVnum5e0
      external RVnum5e0
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

      tcnum1 = RVnum5e0(var,ndim,1.D0,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      RV8e4int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

