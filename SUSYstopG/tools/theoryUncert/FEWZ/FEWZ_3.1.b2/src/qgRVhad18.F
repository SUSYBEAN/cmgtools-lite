CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgRV18e0int(var,ndim,Tcoeff)

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
      double precision qgRVnum14e0
      external qgRVnum14e0
      double precision qgRVnum14e1
      external qgRVnum14e1
      double precision qgRVnum14e2
      external qgRVnum14e2
      double precision qgRVnum14e3
      external qgRVnum14e3
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

      t1 = dlog(z)
      t4 = dlog(1-z)
      t6 = 36*t1
      t7 = 36*t4
      t8 = 1-l1
      t9 = dlog(t8)
      t10 = t9**2
      t12 = t1**2
      t13 = 12*t12
      t14 = 0.3141592653589793D1**2
      t15 = 4*t14
      t16 = t4**2
      t17 = 12*t16
      t19 = t1*t4
      t20 = 24*t19
      t25 = t6-t7+3*t10+t13-t15+t17-18*t9+42-t20+12*t9*t4-12*t9*t1
      t27 = -1/t8
      t33 = dlog(l3)
      t34 = t33**2
      t41 = t6-t7+t13-t15+t17+3*t34-18*t33+42-t20+12*t33*t4-12*t33*t1
      t42 = 1/l3
      tncoff1 = 7.D0/512.D0*t1-7.D0/512.D0*t4-t25*t27/6144+3.D0/512.D0*t
     #12-t14/512+3.D0/512.D0*t16+t41*t42/6144+t12*t1/768-z3/384-t16*t4/7
     #68+15.D0/1024.D0-t4*t12/256-3.D0/256.D0*t19+t16*t1/256-t1*t14/768+
     #t4*t14/768+(-2*t1+2*t4+t9+t33-3)*t27*t42/1024
      t1 = dlog(z)
      t4 = dlog(1-z)
      t6 = dlog(l3)
      t8 = t6**2
      t10 = t4**2
      t12 = t1**2
      t14 = 0.3141592653589793D1**2
      t22 = -36*t1+36*t4+18*t6-3*t8-12*t10-12*t12+4*t14-42+24*t1*t4-12*t
     #6*t4+12*t6*t1
      t23 = 1/l3
      t26 = 1-l1
      t27 = dlog(t26)
      tncoff2 = t22*t23/6144-(-t27+2*t1-2*t4-t6+3)/t26*t23/1024
      t1 = dlog(l1)
      t2 = t1**2
      t5 = dlog(1-z)
      t6 = t5**2
      t8 = 1-l1
      t9 = dlog(t8)
      t10 = t9**2
      t14 = dlog(z)
      t15 = t14**2
      t17 = 0.3141592653589793D1**2
      t33 = -12*t2-12*t6-3*t10+36*t5+18*t9-12*t15+4*t17+36*t1-36*t14-42+
     #12*t9*t14-12*t9*t1-12*t9*t5-24*t5*t1+24*t14*t5+24*t1*t14
      t35 = -1/t8
      t41 = dlog(l3)
      tncoff3 = -t33*t35/6144+(-2*t5-t9-2*t1+2*t14+3-t41)*t35/l3/1024
      t1 = 1-l1
      t2 = dlog(t1)
      t3 = dlog(z)
      t5 = dlog(l1)
      t8 = dlog(1-z)
      t10 = dlog(l3)
      tncoff4 = -(-3+t2-2*t3+2*t5+2*t8+t10)/t1/l3/1024
      t1 = dlog(z)
      t4 = dlog(1-z)
      t6 = 12*t1
      t7 = 12*t4
      t8 = 1-l1
      t9 = dlog(t8)
      t13 = -1/t8
      t16 = dlog(l3)
      t19 = 1/l3
      t22 = t1**2
      t24 = 0.3141592653589793D1**2
      t26 = t4**2
      tncoff5 = 7.D0/1024.D0+3.D0/512.D0*t1-3.D0/512.D0*t4-(18+t6-t7-6*t
     #9)*t13/6144+(18+t6-t7-6*t16)*t19/6144+t22/512-t24/1536+t26/512-t1*
     #t4/256-t13*t19/1024
      t3 = 1/l3
      t6 = dlog(l3)
      t8 = dlog(z)
      t11 = dlog(1-z)
      tncoff6 = 1/(-1+l1)*t3/1024+(-18+6*t6-12*t8+12*t11)*t3/6144
      t1 = 1-l1
      t2 = dlog(t1)
      t4 = dlog(z)
      t6 = dlog(l1)
      t9 = dlog(1-z)
      t13 = -1/t1
      tncoff7 = -(-18+6*t2-12*t4+12*t6+12*t9)*t13/6144+t13/l3/1024
      tncoff8 = -1/(-1+l1)/l3/1024
      t2 = dlog(1-z)
      t4 = dlog(z)
      tncoff9 = 3.D0/1024.D0-t2/512+t4/512+1/l3/1024-1/(-1+l1)/1024
      tncoff10 = -1/l3/1024
      tncoff11 = 1/(-1+l1)/1024
      tncoff12 = 1.D0/1024.D0

CCC Numerator functions

      tcnum1 = qgRVnum14e0(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum14e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum14e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum14e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum14e1(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgRVnum14e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgRVnum14e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgRVnum14e1(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgRVnum14e2(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qgRVnum14e2(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qgRVnum14e2(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qgRVnum14e3(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      qgRV18e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV18e1int(var,ndim,Tcoeff)

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
      double precision qgRVnum14e0
      external qgRVnum14e0
      double precision qgRVnum14e1
      external qgRVnum14e1
      double precision qgRVnum14e2
      external qgRVnum14e2
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
      t6 = dlog(z)
      t9 = 0.3141592653589793D1**2
      t12 = t6**2
      t14 = 1-l1
      t15 = dlog(t14)
      t16 = 2*t2
      t17 = 2*t6
      t20 = -1/t14
      t23 = dlog(l3)
      t25 = 1/l3
      tncoff1 = t3/512-3.D0/512.D0*t2-t6*t2/256+7.D0/1024.D0-t9/1536+3.D
     #0/512.D0*t6+t12/512+(t15+t16-3-t17)*t20/1024-(t23+t16-3-t17)*t25/1
     #024-t20*t25/1024
      t1 = dlog(l3)
      t3 = dlog(1-z)
      t5 = dlog(z)
      t8 = 1/l3
      tncoff2 = -(-t1-2*t3+3+2*t5)*t8/1024+1/(-1+l1)*t8/1024
      t1 = 1-l1
      t2 = dlog(t1)
      t4 = dlog(1-z)
      t6 = dlog(z)
      t8 = dlog(l1)
      t12 = -1/t1
      tncoff3 = (-t2-2*t4+3+2*t6-2*t8)*t12/1024+t12/l3/1024
      tncoff4 = -1/(-1+l1)/l3/1024
      t2 = dlog(1-z)
      t4 = dlog(z)
      tncoff5 = 3.D0/1024.D0-t2/512+t4/512+1/l3/1024-1/(-1+l1)/1024
      tncoff6 = -1/l3/1024
      tncoff7 = 1/(-1+l1)/1024
      tncoff8 = 1.D0/1024.D0

CCC Numerator functions

      tcnum1 = qgRVnum14e0(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum14e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum14e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum14e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgRVnum14e1(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgRVnum14e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgRVnum14e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgRVnum14e2(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff8)

CCC Sum

      qgRV18e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV18e2int(var,ndim,Tcoeff)

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
      double precision qgRVnum14e0
      external qgRVnum14e0
      double precision qgRVnum14e1
      external qgRVnum14e1
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
      t4 = dlog(z)
      tncoff1 = 3.D0/1024.D0-t2/512+t4/512+1/l3/1024-1/(-1+l1)/1024
      tncoff2 = -1/l3/1024
      tncoff3 = 1/(-1+l1)/1024
      tncoff4 = 1.D0/1024.D0

CCC Numerator functions

      tcnum1 = qgRVnum14e0(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgRVnum14e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgRVnum14e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgRVnum14e1(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      qgRV18e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgRV18e3int(var,ndim,Tcoeff)

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
      double precision qgRVnum14e0
      external qgRVnum14e0
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

      tcnum1 = qgRVnum14e0(var,ndim,z,ux,1.D0,l2,0.D0,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      qgRV18e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

