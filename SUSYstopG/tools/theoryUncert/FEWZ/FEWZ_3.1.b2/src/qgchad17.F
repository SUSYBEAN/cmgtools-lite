CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgc17e0int(var,ndim,Tcoeff)

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
      external fillspin,Li2
      double precision qgcnum17e0
      external qgcnum17e0
      double precision qgcnum17e1
      external qgcnum17e1
      double precision qgcnum17e2
      external qgcnum17e2
      double precision qgcnum17e3
      external qgcnum17e3
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
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t1 = 1/l1
      t2 = z*t1
      t4 = dlog(1-z)
      t5 = dlog(z)
      t6 = t4*t5
      t9 = 0.3141592653589793D1**2
      t15 = 1/(-2+l1)
      t19 = 1/l2
      t23 = z*t19
      t26 = t5**2
      t29 = z*t9
      t32 = t4**2
      t37 = dlog(2.D0)
      t38 = t37**2
      t47 = t4*z
      t50 = dlog(l2)
      t51 = t50**2
      t54 = dlog(l1)
      t55 = t54**2
      t58 = -t2*t6/32+z*(2*t9-t9*l1)*t15*t1/1536+t2*t19*t5/32-t23*t6/16+
     #t2*t26/64-t29*t19/768+t23*t32/32+t29*t4/768+t23*t38/128+t2*t32/64+
     #t23*t26/32+t2*t38/256-t47*t26/32+t23*t51/32+t2*t55/256
      t62 = (-2*z+z*l1)*t15
      t63 = t1*t19
      t85 = t23*t5
      t86 = t23*t4
      t92 = t2*t5
      t94 = t2*t4
      t106 = z*t32
      t125 = t62*t63*t37/64-t62*t63*t54/64-t62*t63*t4/32-t62*t63*t50/32+
     #z*t38*t37/768+z*t26*t5/96-z*t32*t4/96+(-t85/16+t86/16)*t50+(-t2*t5
     #4/128+t92/64-t94/64)*t37+(t5*z/128-t47/128)*t38+(z*t26/64-t47*t5/3
     #2+t106/64-t29/1536)*t37+(t106/32-t29/768)*t5-z*z3/384+(-t50*z*t19/
     #32+t85/32-t86/32)*t37+(-t92/64+t94/64)*t54
      tncoff1 = t58+t125
      t4 = l2*z
      t11 = (-z*l1+2*l2+2*z-2*t4-l2*l1+t4*l1)/(-2+l1)
      t12 = 1/l1
      t13 = 1/l2
      t14 = t12*t13
      t16 = dlog(1-z)
      t20 = -l2+t4-z
      t21 = t20*t12
      t23 = dlog(-t20)
      t27 = dlog(2.D0)
      t31 = t20*t13
      t34 = dlog(-z/t20)
      t35 = t34**2
      t38 = t23**2
      t41 = dlog(l2)
      t48 = t41**2
      t51 = t31*t34
      t53 = t31*t16
      t55 = t31*t23
      t65 = 0.3141592653589793D1**2
      t72 = t16**2
      t75 = t27**2
      t78 = dlog(l1)
      tncoff2 = -t11*t14*t16/32+t21*t13*t23/32+t11*t14*t27/64+t31*t35/12
     #8+t31*t38/32-t11*t14*t41/32+t21*t13*t34/64+t31*t48/32+(-t51/32+t53
     #/16-t55/16)*t41+(-t53/32+t55/32)*t34-t31*t23*t16/16+(l2*t65-t4*t65
     #+z*t65)*t13/768+t31*t72/32+t31*t75/128-t11*t14*t78/64+(-t41*t20*t1
     #3/32+t51/64-t53/32+t55/32)*t27
      t1 = -2+l1
      t2 = 1/t1
      t3 = z*t2
      t4 = 1/l1
      t5 = dlog(z)
      t7 = t3*t4*t5
      t10 = dlog(1-z)
      t12 = t3*t4*t10
      t15 = dlog(-t1)
      t16 = t4*t15
      t17 = t3*t16
      t20 = dlog(l1)
      t23 = t4/l2
      t24 = dlog(2.D0)
      t28 = t5**2
      t35 = t10**2
      t39 = dlog(l2)
      t46 = t24**2
      t50 = 0.3141592653589793D1**2
      t55 = t12/16
      t56 = t17/32
      t68 = t15**2
      t78 = t20**2
      tncoff3 = (-t7/32+t12/32+t17/64)*t20+t3*t23*t24/16+t3*t4*t28/32-t3
     #*t23*t20/32+t3*t4*t35/32-t3*t23*t39/16+t3*t23*t5/16+t3*t4*t46/32-z
     #*t50*t2*t4/768+(-t55-t56)*t5+t3*t16*t10/32+(-t3*t4*t20/32+t7/16-t5
     #5-t56)*t24+t3*t4*t68/128-t3*t23*t10/16-t3*t23*t15/32+t3*t4*t78/128
      t2 = -l2+l2*z-z
      t4 = -2+l1
      t5 = 1/t4
      t6 = 2*t2*t5
      t9 = 1/l1/l2
      t10 = dlog(l2)
      t14 = t2*t5
      t17 = dlog(-z/t2)
      t22 = dlog(1-z)
      t27 = dlog(-t4)
      t32 = dlog(-t2)
      t38 = dlog(2.D0)
      t42 = dlog(l1)
      tncoff4 = -t6*t9*t10/32+t14*t9*t17/32-t6*t9*t22/32-t14*t9*t27/32+t
     #14*t9*t32/16+t2*t5*t9*t38/16-t6*t9*t42/64
      t2 = dlog(1-z)
      t3 = z*t2
      t4 = dlog(z)
      t7 = t4**2
      t10 = 1/l1
      t11 = z*t10
      t12 = dlog(l1)
      t15 = dlog(l2)
      t17 = 1/l2
      t20 = z*t17
      t23 = dlog(2.D0)
      t24 = t23**2
      t35 = 0.3141592653589793D1**2
      t44 = t2**2
      tncoff5 = -t3*t4/32+z*t7/64-t11*t12/128-t15*z*t17/32+t20*t4/32+z*t
     #24/256+t11*t4/64-t20*t2/32+(t4*z/64-t3/64)*t23-z*t35/1536-t11*t2/6
     #4+t20*t23/64+t11*t23/128+z*t44/64+(-2*z+z*l1)/(-2+l1)*t10*t17/64
      t1 = dlog(l2)
      t2 = l2*z
      t3 = -l2+t2-z
      t5 = 1/l2
      t8 = t3*t5
      t10 = dlog(-t3)
      t29 = dlog(-z/t3)
      t33 = dlog(1-z)
      t36 = dlog(2.D0)
      tncoff6 = -t1*t3*t5/32+t8*t10/32+(-z*l1+2*l2+2*z-2*t2-l2*l1+t2*l1)
     #/(-2+l1)/l1*t5/64+t8*t29/64-t8*t33/32+t8*t36/64
      t1 = -2+l1
      t3 = z/t1
      t4 = 1/l1
      t6 = dlog(-t1)
      t14 = dlog(l1)
      t18 = dlog(z)
      t23 = dlog(1-z)
      t27 = dlog(2.D0)
      tncoff7 = -t3*t4*t6/64+t3*t4/l2/32-t3*t4*t14/64+t3*t4*t18/32-t3*t4
     #*t23/32+t3*t4*t27/32
      tncoff8 = (-2*l2-2*z+2*l2*z)/(-2+l1)/l1/l2/64
      t2 = dlog(1-z)
      t17 = dlog(2.D0)
      t20 = dlog(z)
      tncoff9 = -z*t2/64+z*(-24+12*l1)/(-2+l1)/l1/1536+z/l2/64+z*t17/128
     #+t20*z/64
      tncoff10 = (-12*z-12*l2+12*l2*z)/l2/768
      tncoff11 = z/(-2+l1)/l1/64
      tncoff12 = z/128

CCC Numerator functions

      tcnum1 = qgcnum17e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgcnum17e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgcnum17e0(var,ndim,z,ux,l1/2,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgcnum17e0(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgcnum17e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgcnum17e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgcnum17e1(var,ndim,z,ux,l1/2,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgcnum17e1(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgcnum17e2(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qgcnum17e2(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qgcnum17e2(var,ndim,z,ux,l1/2,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qgcnum17e3(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      qgc17e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8
     #+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgc17e1int(var,ndim,Tcoeff)

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
      external fillspin,Li2
      double precision qgcnum17e0
      external qgcnum17e0
      double precision qgcnum17e1
      external qgcnum17e1
      double precision qgcnum17e2
      external qgcnum17e2
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
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t2 = dlog(1-z)
      t3 = t2*z
      t4 = dlog(z)
      t7 = t4**2
      t10 = 1/l1
      t11 = z*t10
      t12 = dlog(l1)
      t15 = dlog(l2)
      t17 = 1/l2
      t20 = z*t17
      t23 = dlog(2.D0)
      t24 = t23**2
      t35 = 0.3141592653589793D1**2
      t44 = t2**2
      tncoff1 = -t3*t4/32+z*t7/64-t11*t12/128-t15*z*t17/32+t20*t4/32+z*t
     #24/256+t11*t4/64-t20*t2/32+(t4*z/64-t3/64)*t23-z*t35/1536-t11*t2/6
     #4+t20*t23/64+t11*t23/128+z*t44/64+(-2*z+z*l1)/(-2+l1)*t10*t17/64
      t1 = dlog(l2)
      t2 = l2*z
      t3 = -l2+t2-z
      t5 = 1/l2
      t8 = t3*t5
      t10 = dlog(-t3)
      t29 = dlog(-z/t3)
      t33 = dlog(1-z)
      t36 = dlog(2.D0)
      tncoff2 = -t1*t3*t5/32+t8*t10/32+(-z*l1+2*l2+2*z-2*t2-l2*l1+t2*l1)
     #/(-2+l1)/l1*t5/64+t8*t29/64-t8*t33/32+t8*t36/64
      t1 = -2+l1
      t3 = z/t1
      t4 = 1/l1
      t6 = dlog(-t1)
      t14 = dlog(l1)
      t18 = dlog(z)
      t23 = dlog(1-z)
      t27 = dlog(2.D0)
      tncoff3 = -t3*t4*t6/64+t3*t4/l2/32-t3*t4*t14/64+t3*t4*t18/32-t3*t4
     #*t23/32+t3*t4*t27/32
      tncoff4 = (-2*l2-2*z+2*l2*z)/(-2+l1)/l1/l2/64
      t1 = dlog(2.D0)
      t4 = dlog(z)
      t8 = dlog(1-z)
      tncoff5 = z*t1/128+t4*z/64-t8*z/64+z/l2/64+z/l1/128
      tncoff6 = (-l2+l2*z-z)/l2/64
      tncoff7 = z/(-2+l1)/l1/64
      tncoff8 = z/128

CCC Numerator functions

      tcnum1 = qgcnum17e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgcnum17e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgcnum17e0(var,ndim,z,ux,l1/2,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgcnum17e0(var,ndim,z,ux,l1/2,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgcnum17e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgcnum17e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgcnum17e1(var,ndim,z,ux,l1/2,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgcnum17e2(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff8)

CCC Sum

      qgc17e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgc17e2int(var,ndim,Tcoeff)

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
      external fillspin,Li2
      double precision qgcnum17e0
      external qgcnum17e0
      double precision qgcnum17e1
      external qgcnum17e1
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
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      t1 = dlog(2.D0)
      t4 = dlog(z)
      t8 = dlog(1-z)
      tncoff1 = z*t1/128+t4*z/64-t8*z/64+z/l2/64+z/l1/128
      tncoff2 = (-l2+l2*z-z)/l2/64
      tncoff3 = z/(-2+l1)/l1/64
      tncoff4 = z/128

CCC Numerator functions

      tcnum1 = qgcnum17e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgcnum17e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgcnum17e0(var,ndim,z,ux,l1/2,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgcnum17e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      qgc17e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgc17e3int(var,ndim,Tcoeff)

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
      external fillspin,Li2
      double precision qgcnum17e0
      external qgcnum17e0
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
      call fillspin(1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     & 1.0d0,1.0d0,l1,l2,l3,l4,l5,l6,l7,l8,s1l1,s2l1,s3l1,s4l1,s1l2,
     & s2l2,s3l2,s4l2,c3l1,-999.0d0)
      secsum = 0.0d0


CCC Coefficients

      tncoff1 = z/128

CCC Numerator functions

      tcnum1 = qgcnum17e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      qgc17e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

