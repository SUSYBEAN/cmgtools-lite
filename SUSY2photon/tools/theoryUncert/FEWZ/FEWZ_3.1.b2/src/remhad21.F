CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function rm21e0int(var,ndim,Tcoeff)

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
      double precision rm2num1e0
      external rm2num1e0
      double precision rm2num1e1
      external rm2num1e1
      double precision rm2num1e2
      external rm2num1e2
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

      t1 = k1(z,l1,1.D0,0.D0,l4)
      t2 = -1+z
      t3 = t1*t2
      t4 = 0.3141592653589793D1**2
      t6 = dlog(2.D0)
      t7 = t6**2
      t10 = dlog(1-l4)
      t11 = t10**2
      t13 = dlog(l4)
      t14 = t13**2
      t17 = dlog(1-l1)
      t18 = t17**2
      t21 = dlog(-t2)
      t22 = t21**2
      t24 = dlog(z)
      t25 = t24**2
      t27 = dlog(t1)
      t28 = t27**2
      t30 = dlog(l1)
      t31 = t30**2
      t51 = 2*t4-48*t7-3*t11-3*t14-12*t18-48*t22-12*t25-12*t28-12*t31-12
     #*t17*t13-48*t30*t21-24*t27*t30-48*t27*t21+24*t24*t17-12*t27*t13-12
     #*t27*t10-12*t30*t10+12*t24*t13
      t90 = -24*t6*t13-24*t30*t17-48*t6*t27-96*t6*t21+24*t24*t30+12*t24*
     #t10-24*t21*t10-24*t6*t10-24*t21*t13-24*t27*t17-6*t13*t10-48*t6*t17
     #-12*t30*t13+48*t24*t21-48*t21*t17-12*t17*t10+24*t24*t27+48*t6*t24-
     #48*t6*t30
      t94 = t1*t13
      t95 = t1*t10
      t96 = 1-l2
      t97 = dlog(t96)
      t100 = 2*t1*t17
      t102 = 4*t1*t21
      t104 = 2*t1*t30
      t106 = 2*t1*t27
      t108 = 4*t1*t6
      t110 = 2*t1*t24
      t114 = -1/t96
      t117 = dlog(l3)
      t121 = 1/l3
      tncoff1 = t3*(t51+t90)/3072+t2*(-t94-t95-t97*t1-t100-t102-t104-t10
     #6-t108+t110)*t114/512+t2*(t117*t1+t100+t102+t104+t106+t108-t110+t9
     #5+t94)*t121/512+t3*t114*t121/512
      t1 = -1+z
      t2 = dlog(l3)
      t3 = k1(z,l1,1.D0,l3,l4)
      t6 = dlog(1-l1)
      t10 = dlog(-t1)
      t14 = dlog(1-l3)
      t16 = dlog(2.D0)
      t19 = dlog(l1)
      t22 = dlog(t3)
      t25 = dlog(z)
      t29 = dlog(1-l4)
      t31 = dlog(l4)
      t35 = 1/l3
      tncoff2 = t1*(-t2*t3-2*t6*t3-4*t10*t3-t3*t14-4*t3*t16-2*t3*t19-2*t
     #3*t22+2*t25*t3-t29*t3-t31*t3)*t35/512-t1*t3/(l2-1)*t35/512
      t1 = -1+z
      t2 = dlog(l4)
      t3 = k1(z,l1,l2,0.D0,l4)
      t6 = dlog(1-l4)
      t8 = 1-l2
      t9 = dlog(t8)
      t12 = dlog(1-l1)
      t16 = dlog(-t1)
      t19 = dlog(t3)
      t22 = dlog(2.D0)
      t25 = dlog(z)
      t28 = dlog(l1)
      t31 = dlog(l2)
      t36 = -1/t8
      tncoff3 = t1*(t2*t3+t6*t3+t9*t3+2*t3*t12+4*t3*t16+2*t3*t19+4*t3*t2
     #2-2*t25*t3+2*t3*t28+t31*t3)*t36/512-t1*t3*t36/l3/512
      t2 = k1(z,l1,l2,l3,l4)
      tncoff4 = (-1+z)*t2/(l2-1)/l3/512
      t1 = k1(z,l1,1.D0,0.D0,l4)
      t2 = -1+z
      t3 = t1*t2
      t5 = dlog(1-l4)
      t7 = dlog(z)
      t9 = dlog(t1)
      t11 = dlog(l1)
      t13 = dlog(2.D0)
      t16 = dlog(-t2)
      t19 = dlog(1-l1)
      t21 = dlog(l4)
      tncoff5 = t3*(6*t5-12*t7+12*t9+12*t11+24*t13+24*t16+12*t19+6*t21)/
     #3072+t3/(l2-1)/512-t3/l3/512
      t2 = k1(z,l1,1.D0,l3,l4)
      tncoff6 = (-1+z)*t2/l3/512
      t2 = k1(z,l1,l2,0.D0,l4)
      tncoff7 = -(-1+z)*t2/(l2-1)/512
      t1 = k1(z,l1,1.D0,0.D0,l4)
      tncoff8 = -t1*(-1+z)/512

CCC Numerator functions

      tcnum1 = rm2num1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = rm2num1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = rm2num1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = rm2num1e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,99
     #9.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = rm2num1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = rm2num1e1(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = rm2num1e1(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = rm2num1e2(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff8)

CCC Sum

      rm21e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm21e1int(var,ndim,Tcoeff)

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
      double precision rm2num1e0
      external rm2num1e0
      double precision rm2num1e1
      external rm2num1e1
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

      t1 = k1(z,l1,1.D0,0.D0,l4)
      t2 = -1+z
      t3 = t1*t2
      t5 = dlog(-t2)
      t7 = dlog(t1)
      t10 = dlog(1-l1)
      t12 = dlog(z)
      t14 = dlog(l1)
      t16 = dlog(l4)
      t17 = dlog(2.D0)
      t20 = dlog(1-l4)
      tncoff1 = -t3*(-4*t5-2*t7-2*t10+2*t12-2*t14-t16-4*t17-t20)/512+t3/
     #(l2-1)/512-t3/l3/512
      t2 = k1(z,l1,1.D0,l3,l4)
      tncoff2 = (-1+z)*t2/l3/512
      t2 = k1(z,l1,l2,0.D0,l4)
      tncoff3 = -(-1+z)*t2/(l2-1)/512
      t1 = k1(z,l1,1.D0,0.D0,l4)
      tncoff4 = -t1*(-1+z)/512

CCC Numerator functions

      tcnum1 = rm2num1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = rm2num1e0(var,ndim,z,ux,l1,1.D0,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = rm2num1e0(var,ndim,z,ux,l1,l2,0.D0,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = rm2num1e1(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      rm21e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tn
     #coff4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm21e2int(var,ndim,Tcoeff)

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
      double precision rm2num1e0
      external rm2num1e0
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

      t1 = k1(z,l1,1.D0,0.D0,l4)
      tncoff1 = -t1*(-1+z)/512

CCC Numerator functions

      tcnum1 = rm2num1e0(var,ndim,z,ux,l1,1.D0,0.D0,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      rm21e2int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm21e3int(var,ndim,Tcoeff)

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

      rm21e3int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function rm21e4int(var,ndim,Tcoeff)

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

      rm21e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

