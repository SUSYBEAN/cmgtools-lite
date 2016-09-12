CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgRV3e0int(var,ndim,Tcoeff)

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
      double precision qgRVnum3e0
      external qgRVnum3e0
      double precision qgRVnum3e1
      external qgRVnum3e1
      double precision qgRVnum3e2
      external qgRVnum3e2
      double precision qgRVnum3e3
      external qgRVnum3e3
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
      t2 = dlog(2.D0)
      t6 = dlog(1-z)
      t8 = dlog(z)
      t11 = t6*t2
      t14 = t6*t8
      t16 = t6**2
      t21 = t2*t8
      t26 = t2**2
      t28 = t8**2
      t30 = t1*t2/2048+7.D0/512.D0*t6-7.D0/1024.D0*t8-7.D0/1024.D0*t2+t1
     #1*t8/128+3.D0/256.D0*t14-t16*t8/128+t1*t8/2048-3.D0/512.D0*t21+3.D
     #0/256.D0*t11-t6*t1/1024-3.D0/1024.D0*t26-3.D0/1024.D0*t28
      t42 = dlog(l1)
      t43 = t42*t8
      t46 = t42*t6
      t49 = t42*t2
      t69 = -28+48*t6-24*t8-24*t2-8*t43*l1+16*t46*l1-8*t49*l1-16*t11*l1-
     #16*t14*l1+8*t21*l1+16*t49-12*t42*l1+12*t2*l1+4*t26*l1+12*t8*l1+4*t
     #28*l1
      t75 = t42**2
      t90 = -24*t6*l1+16*t16*l1+16*t43+4*t75*l1-32*t46+32*t14+14*l1+24*t
     #42-l1*t1-16*t21+32*t11-8*t26-8*t28-32*t16-8*t75+2*t1
      t109 = -15.D0/2048.D0-t26*t2/1536+z3/768-t28*t8/1536-3.D0/256.D0*t
     #16+t16*t6/192-(t69+t90)/(-2+l1)/l1/2048+3.D0/4096.D0*t1-t2*t28/512
     #-t26*t8/512+t6*t28/256+t6*t26/256-t16*t2/128
      tncoff1 = t30+t109
      t1 = dlog(2.D0)
      t2 = dlog(z)
      t6 = dlog(1-z)
      t9 = dlog(l1)
      t12 = 2-l1
      t13 = dlog(t12)
      t28 = t2**2
      t30 = 32*t1*t2-32*t6*t2-32*t9*t1+16*t9*t13-32*t1*t13-16*t2*t13-64*
     #t6*t1+32*t9*t6-16*t9*t2+32*t6*t13+8*t28
      t31 = 0.3141592653589793D1**2
      t33 = t6**2
      t35 = t13**2
      t41 = t1**2
      t43 = t9**2
      t46 = 28-2*t31+32*t33+8*t35-48*t6+24*t2-24*t13+48*t1+32*t41+8*t43-
     #24*t9
      tncoff2 = (t30+t46)/t12/l1/2048
      t1 = dlog(l1)
      t6 = dlog(2.D0)
      t9 = dlog(1-z)
      t11 = dlog(z)
      t26 = t11**2
      t28 = 0.3141592653589793D1**2
      t30 = t6**2
      t32 = t9**2
      tncoff3 = -(8*t1+6*l1-4*l1*t1-8*t6+16*t9-8*t11-8*t9*l1+4*t11*l1+4*
     #t6*l1-12)/(-2+l1)/l1/2048-t26/1024+t28/4096-t30/1024-t32/256-3.D0/
     #1024.D0*t6+3.D0/512.D0*t9-3.D0/1024.D0*t11+t9*t6/256+t9*t11/256-t6
     #*t11/512-7.D0/2048.D0
      t2 = dlog(1-z)
      t4 = dlog(z)
      t6 = 2-l1
      t7 = dlog(t6)
      t9 = dlog(2.D0)
      t11 = dlog(l1)
      tncoff4 = (-16*t2+8*t4-8*t7+16*t9-8*t11+12)/t6/l1/2048
      t1 = dlog(z)
      t3 = dlog(2.D0)
      t6 = dlog(1-z)
      tncoff5 = -t1/1024-t3/1024+t6/512-(2*l1-4)/(-2+l1)/l1/2048-3.D0/20
     #48.D0
      tncoff6 = -1/(-2+l1)/l1/512
      tncoff7 = -1.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum3e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum3e0(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum3e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum3e1(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum3e2(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgRVnum3e2(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgRVnum3e3(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff7)

CCC Sum

      qgRV3e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV3e1int(var,ndim,Tcoeff)

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
      double precision qgRVnum3e0
      external qgRVnum3e0
      double precision qgRVnum3e1
      external qgRVnum3e1
      double precision qgRVnum3e2
      external qgRVnum3e2
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
      tncoff1 = -t3/256+t2*t5/256+3.D0/512.D0*t2+t2*t9/256-t12/1024-3.D0
     #/1024.D0*t5-t5*t9/512+t17/4096-7.D0/2048.D0-3.D0/1024.D0*t9-t20/10
     #24+(-4*t22+2*t22*l1+4*t5-2*t5*l1+6-3*l1+4*t9-2*t9*l1-8*t2+4*t2*l1)
     #/(-2+l1)/l1/1024
      t1 = dlog(l1)
      t3 = 2-l1
      t4 = dlog(t3)
      t6 = dlog(2.D0)
      t8 = dlog(z)
      t11 = dlog(1-z)
      tncoff2 = -(4*t1+4*t4-6-8*t6-4*t8+8*t11)/t3/l1/1024
      t2 = dlog(1-z)
      t4 = dlog(2.D0)
      t6 = dlog(z)
      t8 = 2-l1
      tncoff3 = t2/512-t4/1024-3.D0/2048.D0-t6/1024-1/l1/1024
      tncoff4 = -1/(-2+l1)/l1/512
      tncoff5 = -1.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum3e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum3e0(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum3e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum3e1(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum3e2(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff5)

CCC Sum

      qgRV3e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV3e2int(var,ndim,Tcoeff)

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
      double precision qgRVnum3e0
      external qgRVnum3e0
      double precision qgRVnum3e1
      external qgRVnum3e1
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
      tncoff1 = t2/512-t4/1024-3.D0/2048.D0-t6/1024-1/l1/1024
      tncoff2 = -1/(-2+l1)/l1/512
      tncoff3 = -1.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum3e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum3e0(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum3e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff3)

CCC Sum

      qgRV3e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV3e3int(var,ndim,Tcoeff)

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
      double precision qgRVnum3e0
      external qgRVnum3e0
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

      tncoff1 = -1.D0/2048.D0

CCC Numerator functions

      tcnum1 = qgRVnum3e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      qgRV3e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

