CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgRV9e0int(var,ndim,Tcoeff)

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
      double precision qgRVnum7e0
      external qgRVnum7e0
      double precision qgRVnum7e1
      external qgRVnum7e1
      double precision qgRVnum7e2
      external qgRVnum7e2
      double precision qgRVnum7e3
      external qgRVnum7e3
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

      t1 = dlog(l1)
      t3 = dlog(1-z)
      t4 = t1*t3
      t9 = dlog(z)
      t10 = t1*t9
      t14 = t9**2
      t17 = t3**2
      t20 = dlog(2.D0)
      t21 = t20**2
      t28 = t1*t20
      t31 = t1**2
      t34 = 0.3141592653589793D1**2
      t38 = -84+18*t4*l1+36*t1+42*l1+24*t10+36*t9*l1+12*t14*l1+27*t17*l1
     #+3*t21*l1-54*t3*l1-18*t1*l1+12*t28-36*t4+3*t31*l1-t34*l1+18*t20*l1
      t39 = t20*t9
      t42 = t3*t20
      t46 = t3*t9
      t63 = 12*t39*l1-18*t42*l1-6*t31-36*t46*l1-12*t10*l1-6*t28*l1-72*t9
     #-54*t17-36*t20+108*t3+36*t42-24*t39+2*t34-24*t14+72*t46-6*t21
      t87 = 15.D0/1024.D0+(t38+t63)/(-2+l1)/l1/6144-3.D0/512.D0*t42*t9+7
     #.D0/512.D0*t9+t21*t20/6144+27.D0/2048.D0*t17+9.D0/1024.D0*t17*t9+7
     #.D0/1024.D0*t20-21.D0/1024.D0*t3+t14*t9/768-9.D0/1024.D0*t42+3.D0/
     #512.D0*t39-t34*t9/3072
      t109 = t34*t3/2048-3.D0/512.D0*t14*t3+t14*t20/512-t34*t20/6144+9.D
     #0/2048.D0*t17*t20-t34/2048-z3/384+3.D0/512.D0*t14-9.D0/2048.D0*t17
     #*t3+t21*t9/1024-3.D0/2048.D0*t3*t21-9.D0/512.D0*t46+3.D0/2048.D0*t
     #21
      tncoff1 = t87+t109
      t2 = dlog(1-z)
      t3 = 2-l1
      t4 = dlog(t3)
      t7 = dlog(l1)
      t9 = t7**2
      t11 = dlog(z)
      t13 = dlog(2.D0)
      t17 = t11**2
      t19 = t2**2
      t25 = 72*t2*t4-36*t7+6*t9+72*t11+108*t13-72*t4-108*t2+24*t17+54*t1
     #9-72*t13*t4-48*t11*t4
      t34 = t13**2
      t36 = t4**2
      t38 = 0.3141592653589793D1**2
      t46 = 84-72*t2*t11+24*t7*t4+36*t7*t2-36*t7*t13+54*t34+24*t36-2*t38
     #-24*t7*t11-108*t2*t13+72*t13*t11
      tncoff2 = -(t25+t46)/t3/l1/6144
      t2 = dlog(1-z)
      t3 = dlog(2.D0)
      t6 = dlog(z)
      t13 = dlog(l1)
      t35 = t3**2
      t37 = 0.3141592653589793D1**2
      t39 = t6**2
      t41 = t2**2
      tncoff3 = -3.D0/1024.D0*t2*t3-3.D0/512.D0*t2*t6+t3*t6/512+(12*t6*l
     #1+12*t13+6*t3*l1-6*t13*l1-18*t2*l1+18*l1-12*t3-24*t6+36*t2-36)/(-2
     #+l1)/l1/6144+3.D0/1024.D0*t3+3.D0/512.D0*t6-9.D0/1024.D0*t2+t35/20
     #48-t37/6144+t39/512+9.D0/2048.D0*t41+7.D0/1024.D0
      t1 = dlog(l1)
      t3 = dlog(2.D0)
      t5 = 2-l1
      t6 = dlog(t5)
      t9 = dlog(1-z)
      t11 = dlog(z)
      tncoff4 = -(-12*t1+36*t3-24*t6-36*t9+24*t11+36)/t5/l1/6144
      t2 = dlog(1-z)
      t4 = dlog(2.D0)
      t6 = dlog(z)
      tncoff5 = -3.D0/1024.D0*t2+t4/1024+t6/512+3.D0/1024.D0+(-12+6*l1)/
     #(-2+l1)/l1/6144
      tncoff6 = 1/(-2+l1)/l1/512
      tncoff7 = 1.D0/1024.D0

CCC Numerator functions

      tcnum1 = qgRVnum7e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum7e0(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum7e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum7e1(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum7e2(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgRVnum7e2(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgRVnum7e3(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff7)

CCC Sum

      qgRV9e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV9e1int(var,ndim,Tcoeff)

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
      double precision qgRVnum7e0
      external qgRVnum7e0
      double precision qgRVnum7e1
      external qgRVnum7e1
      double precision qgRVnum7e2
      external qgRVnum7e2
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
      t5 = dlog(2.D0)
      t9 = dlog(z)
      t12 = t5**2
      t17 = 0.3141592653589793D1**2
      t20 = t9**2
      t22 = dlog(l1)
      tncoff1 = 9.D0/2048.D0*t3-3.D0/1024.D0*t2*t5-9.D0/1024.D0*t2-3.D0/
     #512.D0*t2*t9+t12/2048+3.D0/1024.D0*t5+t5*t9/512-t17/6144+7.D0/1024
     #.D0+3.D0/512.D0*t9+t20/512-(-2*t22+t22*l1-6*t2+3*t2*l1+2*t5-t5*l1+
     #6-3*l1+4*t9-2*t9*l1)/(-2+l1)/l1/1024
      t1 = dlog(l1)
      t3 = 2-l1
      t4 = dlog(t3)
      t7 = dlog(1-z)
      t9 = dlog(2.D0)
      t11 = dlog(z)
      tncoff2 = (2*t1+4*t4+6*t7-6*t9-6-4*t11)/t3/l1/1024
      t2 = dlog(1-z)
      t4 = dlog(2.D0)
      t6 = dlog(z)
      t8 = 2-l1
      tncoff3 = -3.D0/1024.D0*t2+t4/1024+3.D0/1024.D0+t6/512+1/l1/1024
      tncoff4 = 1/(-2+l1)/l1/512
      tncoff5 = 1.D0/1024.D0

CCC Numerator functions

      tcnum1 = qgRVnum7e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum7e0(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum7e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum7e1(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum7e2(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      qgRV9e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV9e2int(var,ndim,Tcoeff)

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
      double precision qgRVnum7e0
      external qgRVnum7e0
      double precision qgRVnum7e1
      external qgRVnum7e1
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
      t4 = dlog(2.D0)
      t6 = dlog(z)
      tncoff1 = -3.D0/1024.D0*t2+t4/1024+3.D0/1024.D0+t6/512+1/l1/1024
      tncoff2 = 1/(-2+l1)/l1/512
      tncoff3 = 1.D0/1024.D0

CCC Numerator functions

      tcnum1 = qgRVnum7e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum7e0(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum7e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff3)

CCC Sum

      qgRV9e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV9e3int(var,ndim,Tcoeff)

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
      double precision qgRVnum7e0
      external qgRVnum7e0
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

      tncoff1 = 1.D0/1024.D0

CCC Numerator functions

      tcnum1 = qgRVnum7e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      qgRV9e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

