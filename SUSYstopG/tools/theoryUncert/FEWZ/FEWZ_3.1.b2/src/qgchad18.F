CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qgc18e0int(var,ndim,Tcoeff)

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

      t5 = 1/(-2+l1)
      t6 = (-2*z+z*l1)*t5
      t7 = 1/l1
      t8 = 1/l2
      t9 = t7*t8
      t10 = dlog(2.D0)
      t15 = dlog(1-z)
      t19 = z*t7
      t20 = dlog(z)
      t24 = dlog(l1)
      t28 = t20**2
      t31 = t15*z
      t34 = t15**2
      t35 = z*t34
      t37 = 0.3141592653589793D1**2
      t38 = z*t37
      t48 = t10**2
      t52 = t19*t20
      t54 = t19*t15
      t61 = dlog(l2)
      t72 = t15*t20
      t75 = z*t8
      t82 = t6*t9*t10/64-t6*t9*t15/32+t19*t8*t20/32-t6*t9*t24/64+(z*t28/
     #64-t31*t20/32+t35/64-t38/1536)*t10+z*t28*t20/96+(t20*z/128-t31/128
     #)*t48+(-t19*t24/128+t52/64-t54/64)*t10+(-t52/64+t54/64)*t24-t6*t9*
     #t61/32-z*(-2*t37+t37*l1)*t5*t7/1536-t19*t72/32-t75*t72/16+t19*t34/
     #64-t38*t8/768
      t89 = t24**2
      t100 = t61**2
      t118 = t75*t20
      t119 = t75*t15
      t126 = t75*t34/32+t75*t28/32+t19*t48/256+t19*t89/256+t19*t28/64-t3
     #1*t28/32+t38*t15/768+t75*t48/128+t100*z*t8/32-z*t34*t15/96+z*t48*t
     #10/768+(t35/32-t38/768)*t20-z*z3/384+(-t61*z*t8/32+t118/32-t119/32
     #)*t10+(-t118/16+t119/16)*t61
      tncoff1 = t82+t126
      t1 = l2*z
      t2 = -l2+t1-z
      t3 = 1/l2
      t4 = t2*t3
      t6 = dlog(1-z)
      t7 = t4*t6
      t9 = dlog(-t2)
      t10 = t4*t9
      t15 = dlog(-z/t2)
      t17 = t6**2
      t23 = 0.3141592653589793D1**2
      t39 = (-z*l1+2*l2+2*z-2*t1-l2*l1+t1*l1)/(-2+l1)
      t40 = 1/l1
      t41 = t40*t3
      t45 = t2*t40
      t49 = t15**2
      t52 = dlog(l2)
      t56 = t4*t15
      t61 = dlog(2.D0)
      t68 = t61**2
      t71 = t52**2
      t78 = dlog(l1)
      t88 = t9**2
      tncoff2 = (-t7/32+t10/32)*t15+t4*t17/32-t4*t9*t6/16-(-z*t23+t1*t23
     #-l2*t23)*t3/768-t39*t41*t6/32+t45*t3*t9/32+t4*t49/128+(-t52*t2*t3/
     #32+t56/64-t7/32+t10/32)*t61+(-t56/32+t7/16-t10/16)*t52+t4*t68/128+
     #t71*t2*t3/32+t39*t41*t61/64-t39*t41*t78/64-t39*t41*t52/32+t45*t3*t
     #15/64+t4*t88/32
      t1 = -2+l1
      t2 = 1/t1
      t3 = z*t2
      t4 = 1/l1
      t6 = t4/l2
      t7 = dlog(l2)
      t11 = dlog(z)
      t16 = dlog(1-z)
      t21 = dlog(-t1)
      t25 = dlog(2.D0)
      t26 = t25**2
      t30 = dlog(l1)
      t35 = t3*t4*t11
      t38 = t3*t4*t16
      t39 = t38/16
      t40 = t4*t21
      t41 = t40*t3
      t42 = t41/32
      t55 = t30**2
      t59 = 0.3141592653589793D1**2
      t64 = t11**2
      t68 = t16**2
      t72 = t21**2
      tncoff3 = -t3*t6*t7/16+t3*t6*t11/16-t3*t6*t16/16-t3*t6*t21/32+t3*t
     #4*t26/32+(-t3*t4*t30/32+t35/16-t39-t42)*t25+(-t35/32+t38/32+t41/64
     #)*t30+(-t39-t42)*t11+t3*t40*t16/32+t3*t4*t55/128-z*t59*t2*t4/768+t
     #3*t4*t64/32+t3*t4*t68/32+t3*t4*t72/128+t3*t6*t25/16-t3*t6*t30/32
      t2 = -l2+l2*z-z
      t4 = -2+l1
      t5 = 1/t4
      t6 = 2*t2*t5
      t9 = 1/l1/l2
      t11 = dlog(1-z)
      t15 = t2*t5
      t17 = dlog(-t4)
      t22 = dlog(-t2)
      t28 = dlog(2.D0)
      t32 = dlog(l1)
      t36 = dlog(l2)
      t42 = dlog(-z/t2)
      tncoff4 = -t6*t9*t11/32-t15*t9*t17/32+t15*t9*t22/16+t2*t5*t9*t28/1
     #6-t6*t9*t32/64-t6*t9*t36/32+t15*t9*t42/32
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
      t5 = dlog(1-z)
      t8 = dlog(2.D0)
      t11 = dlog(z)
      tncoff9 = z/l2/64-t5*z/64+z*t8/128+t11*z/64-z*(24-12*l1)/(-2+l1)/l
     #1/1536
      tncoff10 = -(12*z+12*l2-12*l2*z)/l2/768
      tncoff11 = z/(-2+l1)/l1/64
      tncoff12 = z/128

CCC Numerator functions

      tcnum1 = qgcnum17e0(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgcnum17e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgcnum17e0(var,ndim,z,ux,1-l1/2,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgcnum17e0(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgcnum17e1(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgcnum17e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgcnum17e1(var,ndim,z,ux,1-l1/2,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgcnum17e1(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qgcnum17e2(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qgcnum17e2(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qgcnum17e2(var,ndim,z,ux,1-l1/2,0.D0,l3,l4,l5,l6,l7,l8,9
     #99.D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qgcnum17e3(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff12)

CCC Sum

      qgc18e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8
     #+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgc18e1int(var,ndim,Tcoeff)

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
      tncoff5 = z*t1/128+t4*z/64-z*t8/64+z/l2/64+z/l1/128
      tncoff6 = (-l2+l2*z-z)/l2/64
      tncoff7 = z/(-2+l1)/l1/64
      tncoff8 = z/128

CCC Numerator functions

      tcnum1 = qgcnum17e0(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgcnum17e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgcnum17e0(var,ndim,z,ux,1-l1/2,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgcnum17e0(var,ndim,z,ux,1-l1/2,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qgcnum17e1(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qgcnum17e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qgcnum17e1(var,ndim,z,ux,1-l1/2,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qgcnum17e2(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff8)

CCC Sum

      qgc18e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgc18e2int(var,ndim,Tcoeff)

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

      tcnum1 = qgcnum17e0(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qgcnum17e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qgcnum17e0(var,ndim,z,ux,1-l1/2,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qgcnum17e1(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff4)

CCC Sum

      qgc18e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*t
     #ncoff4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qgc18e3int(var,ndim,Tcoeff)

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

      tcnum1 = qgcnum17e0(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff1)

CCC Sum

      qgc18e3int = tcnum1*tncoff1

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

