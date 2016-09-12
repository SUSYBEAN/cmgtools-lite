CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function ggrm12e0int(var,ndim,Tcoeff)

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
      double precision ggrmnum4e0
      external ggrmnum4e0
      double precision ggrmnum4e1
      external ggrmnum4e1
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
      t2 = t1**2
      t3 = l3-1
      t4 = t2*t3
      t6 = dlog(1-l1)
      t7 = dlog(z)
      t10 = dlog(1-l4)
      t11 = dlog(l4)
      t12 = dlog(l3)
      t15 = dlog(-t3)
      t17 = f2(l3,z)
      t18 = dlog(t17)
      t19 = f3(l3,z)
      t22 = dlog(1-l1*t19)
      t23 = dlog(2.D0)
      t26 = dlog(-t1)
      t28 = dlog(l1)
      t29 = t6-2-2*t7+t10+t11+2*t12+2*t15-t18+t22+4*t23+4*t26+t28
      t30 = s23r2(z,l1,0.D0,l3,l4)
      t31 = 1/t30
      t33 = f1(l3,z,l1)
      tncoff1 = 9.D0/32768.D0*t4*t29*t31*t33-9.D0/32768.D0*t4/l8*t31*t33
     #-9.D0/32768.D0*t4/l2*t31*t33
      t1 = -1+z
      t2 = t1**2
      t3 = l3-1
      t4 = t2*t3
      t6 = dlog(1-l1)
      t7 = dlog(z)
      t10 = dlog(1-l4)
      t11 = dlog(l4)
      t12 = dlog(l3)
      t15 = dlog(-t3)
      t17 = f2(l3,z)
      t18 = dlog(t17)
      t19 = f3(l3,z)
      t22 = dlog(1-l1*t19)
      t23 = dlog(2.D0)
      t26 = dlog(-t1)
      t28 = dlog(l1)
      t29 = t6-2-2*t7+t10+t11+2*t12+2*t15-t18+t22+4*t23+4*t26+t28
      t30 = s23r2(z,l1,0.D0,l3,l4)
      t31 = 1/t30
      t33 = f1(l3,z,l1)
      tncoff2 = 9.D0/32768.D0*t4*t29*t31*t33+9.D0/32768.D0*t4/(-1+l8)*t3
     #1*t33-9.D0/32768.D0*t4/l2*t31*t33
      t2 = (-1+z)**2
      t9 = s23r2(z,l1,0.D0,l3,l4)
      t12 = f1(l3,z,l1)
      tncoff3 = -9.D0/32768.D0*t2*(l3-1)/l8/(-1+l8)/t9*t12
      t3 = (-1+z)**2
      t6 = s23r2(z,l1,l2,l3,l4)
      t9 = f1(l3,z,l1)
      tncoff4 = 9.D0/32768.D0*(l3-1)*t3/l2/t6*t9
      t3 = (-1+z)**2
      t6 = s23r2(z,l1,l2,l3,l4)
      t9 = f1(l3,z,l1)
      tncoff5 = 9.D0/32768.D0*(l3-1)*t3/l2/t6*t9
      t2 = (-1+z)**2
      t5 = s23r2(z,l1,0.D0,l3,l4)
      t7 = f1(l3,z,l1)
      tncoff6 = -9.D0/32768.D0*t2*(l3-1)/t5*t7
      t2 = (-1+z)**2
      t5 = s23r2(z,l1,0.D0,l3,l4)
      t7 = f1(l3,z,l1)
      tncoff7 = -9.D0/32768.D0*t2*(l3-1)/t5*t7

CCC Numerator functions

      tcnum1 = ggrmnum4e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ggrmnum4e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = ggrmnum4e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = ggrmnum4e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,0.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = ggrmnum4e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,1.D0,-999.D
     #0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = ggrmnum4e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = ggrmnum4e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff7)

CCC Sum

      ggrm12e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function ggrm12e1int(var,ndim,Tcoeff)

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
      double precision ggrmnum4e0
      external ggrmnum4e0
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

      t2 = (-1+z)**2
      t5 = s23r2(z,l1,0.D0,l3,l4)
      t7 = f1(l3,z,l1)
      tncoff1 = -9.D0/32768.D0*t2*(l3-1)/t5*t7
      t2 = (-1+z)**2
      t5 = s23r2(z,l1,0.D0,l3,l4)
      t7 = f1(l3,z,l1)
      tncoff2 = -9.D0/32768.D0*t2*(l3-1)/t5*t7

CCC Numerator functions

      tcnum1 = ggrmnum4e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,0.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = ggrmnum4e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,1.D0,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      ggrm12e1int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function ggrm12e2int(var,ndim,Tcoeff)

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

      ggrm12e2int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

