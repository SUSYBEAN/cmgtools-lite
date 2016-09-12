CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgsc22e0int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,Li2,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      double precision mur,muf,Lf,m2,scm
      COMMON/scale/ mur,muf
      COMMON/params/ m2,scm
      external fillspin,Li2
      double precision qgscnum21e0
      external qgscnum21e0
      double precision qgscnum21e1
      external qgscnum21e1
      double precision qgscnum21e2
      external qgscnum21e2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      Lf=dlog(muf**2/m2)
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t1 = Lf*z
      t2 = dlog(2.D0)
      t3 = t2**2
      t6 = dlog(z)
      t9 = dlog(1-z)
      t14 = t6**2
      t20 = t9**2
      t23 = 0.3141592653589793D1**2
      t26 = 1/l2
      t30 = dlog(l2)
      t40 = 1/l1
      t44 = dlog(l1)
      tncoff1 = -t1*t3/256+(-t1*t6/64+t1*t9/64)*t2-t1*t14/64+t1*t9*t6/32
     #-t1*t20/64+t1*t23/1536-t1*t26*t2/64+t1*t26*t30/32-t1*t26*t6/32+t1*
     #t26*t9/32-t1*t40*t2/128+t44*z*Lf*t40/128-t1*t40*t6/64+t1*t40*t9/64
     #-Lf*(-2*z+z*l1)/(-2+l1)*t40*t26/64
      t1 = l2*z
      t2 = -l2+t1-z
      t3 = t2*Lf
      t4 = 1/l2
      t5 = dlog(2.D0)
      t9 = dlog(l2)
      t15 = dlog(-z/t2)
      t20 = dlog(1-z)
      t25 = dlog(-t2)
      tncoff2 = -t3*t4*t5/64+t3*t4*t9/32-t3*t4*t15/64+t3*t4*t20/32-t3*t4
     #*t25/32-Lf*(2*l2-l2*l1-2*t1+t1*l1+2*z-z*l1)/(-2+l1)/l1*t4/64
      t1 = Lf*z
      t2 = -2+l1
      t3 = 1/t2
      t4 = 1/l1
      t5 = t3*t4
      t6 = dlog(2.D0)
      t10 = dlog(l1)
      t16 = dlog(z)
      t21 = dlog(1-z)
      t26 = dlog(-t2)
      tncoff3 = -t1*t5*t6/32+t10*z*Lf*t3*t4/64-t1*t5*t16/32+t1*t5*t21/32
     #+t5*t1*t26/64-t1*t5/l2/32
      tncoff4 = -Lf*(-2*l2+2*l2*z-2*z)/(-2+l1)/l1/l2/64
      t1 = Lf*z
      t2 = dlog(2.D0)
      t5 = dlog(z)
      t9 = dlog(1-z)
      tncoff5 = -t1*t2/128-t1*t5/64+t1*t9/64-t1/l2/64-t1/l1/128
      tncoff6 = -Lf*(-l2+l2*z-z)/l2/64
      tncoff7 = -Lf*z/(-2+l1)/l1/64
      tncoff8 = -Lf*z/128

CCC Numerator functions

      tcnum1 = qgscnum21e0(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgscnum21e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgscnum21e0(var,ndim,z,ux,1-l1/2,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgscnum21e0(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgscnum21e1(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgscnum21e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgscnum21e1(var,ndim,z,ux,1-l1/2,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgscnum21e2(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff8)

CCC Sum

      qgsc22e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgsc22e1int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,Li2,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      double precision mur,muf,Lf,m2,scm
      COMMON/scale/ mur,muf
      COMMON/params/ m2,scm
      external fillspin,Li2
      double precision qgscnum21e0
      external qgscnum21e0
      double precision qgscnum21e1
      external qgscnum21e1
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      Lf=dlog(muf**2/m2)
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t1 = Lf*z
      t2 = dlog(2.D0)
      t5 = dlog(z)
      t9 = dlog(1-z)
      tncoff1 = -t1*t2/128-t1*t5/64+t1*t9/64-t1/l2/64-t1/l1/128
      tncoff2 = -Lf*(-l2+l2*z-z)/l2/64
      tncoff3 = -Lf*z/(-2+l1)/l1/64
      tncoff4 = -Lf*z/128

CCC Numerator functions

      tcnum1 = qgscnum21e0(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgscnum21e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgscnum21e0(var,ndim,z,ux,1-l1/2,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgscnum21e1(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      qgsc22e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgsc22e2int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,Li2,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      double precision mur,muf,Lf,m2,scm
      COMMON/scale/ mur,muf
      COMMON/params/ m2,scm
      external fillspin,Li2
      double precision qgscnum21e0
      external qgscnum21e0
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      Lf=dlog(muf**2/m2)
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      tncoff1 = -Lf*z/128

CCC Numerator functions

      tcnum1 = qgscnum21e0(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      qgsc22e2int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgsc22e3int(var,ndim,Tcoeff)

      IMPLICIT double precision(s,t)
      double precision Tcoeff
      INTEGER ndim
      double precision var(ndim)
      double precision l1,l2,l3,l4,l5,l6,l7,l8,z,u,t,r,ux,Pi,buf,jac
      double precision secsum,z3,NF,Li2,e1temp
      double precision k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,
     & s34r2,s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      double precision Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322,Fjet
      external k1,s13f,s23f,s14f,s24f,s3vf,s4vf,f1,f2,f3,s34r2,
     & s23r2,s1vr2,s2vr2,s23rf,s34rf,s1vrf,s13rf
      external Fjet,Diquad71,Diquad721,Diquad722,
     & Dsquar311,Dsquar312,Dsquar321,Dsquar322
      double precision s1l1,s2l1,s3l1,s4l1,s1l2,s2l2,s3l2,s4l2,c3l1
      double precision yp,zp
      double precision mur,muf,Lf,m2,scm
      COMMON/scale/ mur,muf
      COMMON/params/ m2,scm
      external fillspin,Li2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l5 = jac*var(7)+buf
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = jac*var(9)+buf
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      yp = l2*(1.0d0-z)+z
      zp=z/yp
      Lf=dlog(muf**2/m2)
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients


CCC Numerator functions


CCC Sum

      qgsc22e3int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

