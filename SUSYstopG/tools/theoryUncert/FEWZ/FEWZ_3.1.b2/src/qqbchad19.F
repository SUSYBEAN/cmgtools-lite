CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function qqbc19e0int(var,ndim,Tcoeff)

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
      double precision qqbcnum19e0
      external qqbcnum19e0
      double precision qqbcnum19e1
      external qqbcnum19e1
      double precision qqbcnum19e2
      external qqbcnum19e2
      double precision qqbcnum19e3
      external qqbcnum19e3
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

      t1 = 0.3141592653589793D1**2
      t2 = z*t1
      t3 = 1/l2
      t6 = z*t3
      t7 = dlog(z)
      t8 = t7**2
      t15 = 1/(-1+l1)
      t16 = 1/l1
      t17 = t15*t16
      t20 = z*t16
      t21 = dlog(l1)
      t22 = t21**2
      t27 = dlog(1-z)
      t32 = t27*t7
      t37 = t27**2
      t48 = -t2*t3/768+t6*t8/32-z*(-t1+t1*l1)*t17/1536+t20*t22/256+(-t20
     #*t7/64+t20*t27/64)*t21-t20*t32/32+t20*t8/64+t6*t37/32+t20*t37/64+z
     #*t8*t7/96-z*t37*t27/96
      t60 = -z+z*l1
      t61 = t60*t15
      t62 = t16*t3
      t68 = dlog(l2)
      t73 = t68**2
      t90 = -t27*z*t8/32+(z*t37/32-t2/768)*t7-z*z3/384-t61*t62*t21/64+t2
     #*t27/768-t68*t60*t17*t3/32+t73*z*t3/32+t20*t3*t7/32-t61*t62*t27/32
     #+(-t6*t7/16+t6*t27/16)*t68-t6*t32/16
      tncoff1 = t48+t90
      t1 = 0.3141592653589793D1**2
      t2 = t1*l2
      t6 = 1/l2
      t9 = l2*z
      t10 = -l2+t9-z
      t11 = t10*t6
      t13 = dlog(-t10)
      t15 = dlog(1-z)
      t19 = t15**2
      t24 = dlog(-z/t10)
      t25 = t24**2
      t28 = t13**2
      t31 = t11*t15
      t32 = t11*t13
      t39 = -t9+z-z*l1+l2+t9*l1-l2*l1
      t41 = 1/(-1+l1)
      t42 = t39*t41
      t43 = 1/l1
      t44 = t43*t6
      t45 = dlog(l1)
      t49 = dlog(l2)
      t55 = t49**2
      t59 = t10*t43
      tncoff2 = -(t2*z-z*t1-t2)*t6/768-t11*t13*t15/16+t11*t19/32+t11*t25
     #/128+t11*t28/32+(-t31/32+t32/32)*t24-t42*t44*t45/64-t49*t39*t41*t4
     #3*t6/32+t55*t10*t6/32+t59*t6*t24/64+(-t11*t24/32+t31/16-t32/16)*t4
     #9-t42*t44*t15/32+t59*t6*t13/32
      t1 = -1+l1
      t2 = 1/t1
      t3 = t2*z
      t5 = dlog(-t1)
      t9 = dlog(1-z)
      t13 = dlog(z)
      t15 = t13**2
      t21 = 1/l2
      t22 = z*t21
      t27 = dlog(l2)
      t32 = t5**2
      t35 = t9**2
      t40 = 0.3141592653589793D1**2
      t41 = z*t40
      t50 = (t3*t5/64+t3*t9/32)*t13-t3*t15/64-t3*t9*t5/64+(-t22*t13/16+t
     #22*t9/16)*t27-t22*t9*t13/16-t3*t32/256-t3*t35/64+t22*t35/32-t41*t2
     #1/768+t27*z*t2*t21/32+t41*t2/1536
      t80 = t27**2
      t84 = t22*t15/32-t3*t21*t13/32+t3*t21*t5/64+t3*t21*t9/32+t41*t9/76
     #8+z*t15*t13/96-z*t35*t9/96-t9*z*t15/32+(z*t35/32-t41/768)*t13-z*z3
     #/384+t80*z*t21/32
      tncoff3 = t50+t84
      t1 = l2*z
      t2 = -l2+t1-z
      t3 = 1/l2
      t4 = t2*t3
      t7 = dlog(-z/t2)
      t11 = dlog(1-z)
      t12 = t4*t11
      t15 = dlog(-t2)
      t16 = t4*t15
      t19 = dlog(l2)
      t27 = 0.3141592653589793D1**2
      t28 = t27*l2
      t34 = t7**2
      t37 = t11**2
      t43 = l2*l1-t1*l1+z*l1
      t45 = -1+l1
      t46 = 1/t45
      t47 = 1/l1
      t52 = t2*t46
      t59 = t43*t46
      t60 = t47*t3
      t62 = dlog(-t45)
      t69 = t15**2
      t72 = t19**2
      tncoff4 = (-t4*t7/32+t12/16-t16/16)*t19+(-t12/32+t16/32)*t7-t4*t15
     #*t11/16-(t28*z-z*t27-t28)*t3/768+t4*t34/128+t4*t37/32-t19*t43*t46*
     #t47*t3/32-t52*t3*t15/32-t52*t3*t7/64-t59*t60*t62/64-t59*t60*t11/32
     #+t4*t69/32+t72*t2*t3/32
      t1 = -1+l1
      t2 = 1/t1
      t3 = t2*z
      t4 = 1/l1
      t5 = 1/l2
      t6 = t4*t5
      t7 = dlog(l1)
      t11 = dlog(l2)
      t13 = t2*t4
      t17 = dlog(z)
      t22 = dlog(-t1)
      t27 = dlog(1-z)
      t31 = t22**2
      t39 = t3*t4*t22
      t41 = t4*t27
      t42 = t3*t41
      t46 = t17**2
      t57 = t7**2
      t61 = t27**2
      t65 = 0.3141592653589793D1**2
      tncoff5 = -t3*t6*t7/64-t11*z*t13*t5/32+t3*t6*t17/32-t3*t6*t22/64-t
     #3*t6*t27/32+t3*t4*t31/256+(-t3*t4*t17/64+t39/128+t42/64)*t7+t3*t4*
     #t46/64+(-t39/64-t42/32)*t17+t3*t41*t22/64+t3*t4*t57/256+t3*t4*t61/
     #64-z*t65*t13/1536
      t2 = -l2+l2*z-z
      t3 = -1+l1
      t4 = 1/t3
      t5 = t2*t4
      t6 = 1/l1
      t7 = 1/l2
      t8 = t6*t7
      t9 = dlog(l1)
      t13 = dlog(l2)
      t21 = dlog(-z/t2)
      t26 = dlog(-t3)
      t31 = dlog(1-z)
      t36 = dlog(-t2)
      tncoff6 = -t5*t8*t9/64-t13*t2*t4*t6*t7/32+t5*t8*t21/64-t5*t8*t26/6
     #4-t5*t8*t31/32+t5*t8*t36/32
      t1 = 1/l2
      t2 = z*t1
      t4 = dlog(1-z)
      t7 = 1/l1
      t8 = z*t7
      t9 = dlog(l1)
      t12 = dlog(z)
      t15 = t12**2
      t23 = t4**2
      t34 = 0.3141592653589793D1**2
      t37 = dlog(l2)
      tncoff7 = -t2*t4/32-t8*t9/128+t8*t12/64+z*t15/64-t4*z*t12/32-t8*t4
     #/64+z*t23/64+(-z+z*l1)/(-1+l1)*t7*t1/64-z*t34/1536-t2*t37/32+t2*t1
     #2/32
      t1 = l2*z
      t2 = -l2+t1-z
      t3 = 1/l2
      t4 = t2*t3
      t6 = dlog(1-z)
      t11 = dlog(-z/t2)
      t15 = dlog(-t2)
      t18 = dlog(l2)
      tncoff8 = -t4*t6/32+t4*t11/64+t4*t15/32-t4*t18/32+(-t1+z-z*l1+l2+t
     #1*l1-l2*l1)/(-1+l1)/l1*t3/64
      t2 = dlog(1-z)
      t3 = t2**2
      t6 = -1+l1
      t8 = z/t6
      t9 = dlog(z)
      t13 = dlog(-t6)
      t16 = 1/l2
      t17 = z*t16
      t18 = dlog(l2)
      t29 = t9**2
      t35 = 0.3141592653589793D1**2
      tncoff9 = z*t3/64-t8*t9/64+t8*t13/128-t17*t18/32+t17*t9/32-t17*t2/
     #32+t8*t2/64-t8*t16/64+z*t29/64-t2*z*t9/32-z*t35/1536
      t1 = l2*z
      t2 = -l2+t1-z
      t3 = 1/l2
      t4 = t2*t3
      t6 = dlog(-t2)
      t20 = dlog(l2)
      t25 = dlog(-z/t2)
      t29 = dlog(1-z)
      tncoff10 = t4*t6/32+(l2*l1-t1*l1+z*l1)/(-1+l1)/l1*t3/64-t4*t20/32+
     #t4*t25/64-t4*t29/32
      t1 = -1+l1
      t3 = z/t1
      t4 = 1/l1
      t5 = dlog(l1)
      t9 = dlog(z)
      t14 = dlog(-t1)
      t19 = dlog(1-z)
      tncoff11 = -t3*t4*t5/128+t3*t4*t9/64-t3*t4*t14/128-t3*t4*t19/64+t3
     #*t4/l2/64
      tncoff12 = (-l2+l2*z-z)/(-1+l1)/l1/l2/64
      t1 = dlog(z)
      t5 = dlog(1-z)
      t8 = 1-l1
      tncoff13 = z*t1/64-t5*z/64+z/l1/128+z/l2/64
      tncoff14 = -(12*l2+12*z-12*l2*z)/l2/768
      t1 = dlog(z)
      t5 = dlog(1-z)
      tncoff15 = z*t1/64-t5*z/64+z/l2/64-z/(-1+l1)/128
      tncoff16 = -(12*l2+12*z-12*l2*z)/l2/768
      tncoff17 = z/(-1+l1)/l1/128
      tncoff18 = z/128
      tncoff19 = z/128

CCC Numerator functions

      tcnum1 = qqbcnum19e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum19e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbcnum19e0(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbcnum19e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbcnum19e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbcnum19e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbcnum19e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qqbcnum19e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qqbcnum19e1(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qqbcnum19e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qqbcnum19e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qqbcnum19e1(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0
     #,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qqbcnum19e2(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = qqbcnum19e2(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = qqbcnum19e2(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = qqbcnum19e2(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = qqbcnum19e2(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = qqbcnum19e3(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = qqbcnum19e3(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff19)

CCC Sum

      qqbc19e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncof
     #f16+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc19e1int(var,ndim,Tcoeff)

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
      double precision qqbcnum19e0
      external qqbcnum19e0
      double precision qqbcnum19e1
      external qqbcnum19e1
      double precision qqbcnum19e2
      external qqbcnum19e2
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

      t1 = 1/l2
      t2 = z*t1
      t4 = dlog(1-z)
      t7 = 1/l1
      t8 = z*t7
      t9 = dlog(l1)
      t12 = dlog(z)
      t15 = t12**2
      t23 = t4**2
      t34 = 0.3141592653589793D1**2
      t37 = dlog(l2)
      tncoff1 = -t2*t4/32-t8*t9/128+t8*t12/64+z*t15/64-t4*z*t12/32-t8*t4
     #/64+z*t23/64+(-z+z*l1)/(-1+l1)*t7*t1/64-z*t34/1536-t2*t37/32+t2*t1
     #2/32
      t1 = l2*z
      t2 = -l2+t1-z
      t3 = 1/l2
      t4 = t2*t3
      t6 = dlog(1-z)
      t11 = dlog(-z/t2)
      t15 = dlog(-t2)
      t18 = dlog(l2)
      tncoff2 = -t4*t6/32+t4*t11/64+t4*t15/32-t4*t18/32+(-t1+z-z*l1+l2+t
     #1*l1-l2*l1)/(-1+l1)/l1*t3/64
      t2 = dlog(1-z)
      t3 = t2**2
      t6 = -1+l1
      t8 = z/t6
      t9 = dlog(z)
      t13 = dlog(-t6)
      t16 = 1/l2
      t17 = z*t16
      t18 = dlog(l2)
      t29 = t9**2
      t35 = 0.3141592653589793D1**2
      tncoff3 = z*t3/64-t8*t9/64+t8*t13/128-t17*t18/32+t17*t9/32-t17*t2/
     #32+t8*t2/64-t8*t16/64+z*t29/64-t2*z*t9/32-z*t35/1536
      t1 = l2*z
      t2 = -l2+t1-z
      t3 = 1/l2
      t4 = t2*t3
      t6 = dlog(-t2)
      t20 = dlog(l2)
      t25 = dlog(-z/t2)
      t29 = dlog(1-z)
      tncoff4 = t4*t6/32+(l2*l1-t1*l1+z*l1)/(-1+l1)/l1*t3/64-t4*t20/32+t
     #4*t25/64-t4*t29/32
      t1 = -1+l1
      t3 = z/t1
      t4 = 1/l1
      t5 = dlog(l1)
      t9 = dlog(z)
      t14 = dlog(-t1)
      t19 = dlog(1-z)
      tncoff5 = -t3*t4*t5/128+t3*t4*t9/64-t3*t4*t14/128-t3*t4*t19/64+t3*
     #t4/l2/64
      tncoff6 = (-l2+l2*z-z)/(-1+l1)/l1/l2/64
      t1 = dlog(z)
      t8 = dlog(1-z)
      tncoff7 = z*t1/64+z/l2/64-t8*z/64+z/l1/128
      tncoff8 = (-l2+l2*z-z)/l2/64
      t1 = dlog(z)
      t5 = dlog(1-z)
      tncoff9 = z*t1/64-t5*z/64+z/l2/64-z/(-1+l1)/128
      tncoff10 = (-l2+l2*z-z)/l2/64
      tncoff11 = z/(-1+l1)/l1/128
      tncoff12 = z/128
      tncoff13 = z/128

CCC Numerator functions

      tcnum1 = qqbcnum19e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum19e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbcnum19e0(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbcnum19e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbcnum19e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbcnum19e0(var,ndim,z,ux,l1,l2,l3,l4,l5,l6,l7,l8,999.D0,
     #999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbcnum19e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = qqbcnum19e1(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = qqbcnum19e1(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = qqbcnum19e1(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = qqbcnum19e1(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,999.
     #D0,999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = qqbcnum19e2(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = qqbcnum19e2(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,99
     #9.D0,999.D0,-999.D0,Tcoeff*tncoff13)

CCC Sum

      qqbc19e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff
     #8+tcnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff1
     #2+tcnum13*tncoff13

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc19e2int(var,ndim,Tcoeff)

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
      double precision qqbcnum19e0
      external qqbcnum19e0
      double precision qqbcnum19e1
      external qqbcnum19e1
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

      t1 = dlog(z)
      t8 = dlog(1-z)
      tncoff1 = z*t1/64+z/l2/64-t8*z/64+z/l1/128
      tncoff2 = (-l2+l2*z-z)/l2/64
      t1 = dlog(z)
      t5 = dlog(1-z)
      tncoff3 = z*t1/64-t5*z/64+z/l2/64-z/(-1+l1)/128
      tncoff4 = (-l2+l2*z-z)/l2/64
      tncoff5 = z/(-1+l1)/l1/128
      tncoff6 = z/128
      tncoff7 = z/128

CCC Numerator functions

      tcnum1 = qqbcnum19e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum19e0(var,ndim,z,ux,0.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = qqbcnum19e0(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = qqbcnum19e0(var,ndim,z,ux,1.D0,l2,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = qqbcnum19e0(var,ndim,z,ux,l1,0.D0,l3,l4,l5,l6,l7,l8,999.D
     #0,999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = qqbcnum19e1(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = qqbcnum19e1(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff7)

CCC Sum

      qqbc19e2int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*
     #tncoff4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc19e3int(var,ndim,Tcoeff)

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
      double precision qqbcnum19e0
      external qqbcnum19e0
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
      tncoff2 = z/128

CCC Numerator functions

      tcnum1 = qqbcnum19e0(var,ndim,z,ux,0.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = qqbcnum19e0(var,ndim,z,ux,1.D0,0.D0,l3,l4,l5,l6,l7,l8,999
     #.D0,999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      qqbc19e3int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function qqbc19e4int(var,ndim,Tcoeff)

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


CCC Numerator functions


CCC Sum

      qqbc19e4int = 0

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

