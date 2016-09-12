CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function CD2e0int(var,ndim,Tcoeff)

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
      double precision CDnum1e0
      external CDnum1e0
      double precision CDnum1e1
      external CDnum1e1
      double precision CDnum1e2
      external CDnum1e2
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l5 = var(7)
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = var(9)
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t1 = -1+z
      t2 = 1/l8
      t3 = t1*t2
      t5 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t6 = 1/t5
      t7 = 1/l3*t6
      t9 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t10 = 1/t9
      t11 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t12 = 1/t11
      t14 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t15 = 1/t14
      t16 = f1(1.D0,z,1.D0)
      t18 = t10*t12*t15*t16
      t21 = 1/l2
      t22 = t21*t6
      t27 = dlog(1-l4)
      t28 = t27**2
      t31 = dlog(-t1)
      t32 = t31**2
      t34 = dlog(l1)
      t35 = t34**2
      t37 = dlog(l4)
      t38 = t37**2
      t40 = f3(1.D0,z)
      t42 = dlog(1-t40)
      t49 = dlog(z)
      t52 = dlog(2.D0)
      t57 = t52**2
      t59 = 0.3141592653589793D1**2
      t60 = f2(1.D0,z)
      t61 = dlog(t60)
      t62 = t61**2
      t66 = t42**2
      t74 = 3*t28+48*t32+3*t35+3*t38+6*t27*t42+6*t34*t42+6*t37*t42-48*t3
     #1*t49+48*t31*t52+24*t37*t31+12*t57-t59+3*t62-24*t31*t61+3*t66-12*t
     #37*t49-12*t34*t49+24*t34*t31
      t77 = t49**2
      t113 = 6*t27*t37+12*t77+12*t42*t52+24*t31*t42-12*t42*t49-6*t42*t61
     #+12*t61*t49-6*t34*t61-12*t61*t52-6*t37*t61+12*t37*t52+6*t37*t34-12
     #*t27*t49+24*t27*t31+12*t27*t52+6*t27*t34+12*t34*t52-24*t49*t52-6*t
     #27*t61
      t119 = 2*t52
      t120 = 4*t31
      t121 = 2*t49
      t122 = dlog(l2)
      t129 = dlog(l3)
      t143 = dlog(l8)
      t149 = t27+t37+t34+t119+t120+t42-t61-t121-4*t31*l8-t34*l8-t27*l8+2
     #*t49*l8-t42*l8+t143-t37*l8+t61*l8-2*t52*l8-t143*l8
      tncoff1 = -t3*t7*t18/1024-t3*t22*t18/1024-t1*(t74+t113)*t6*t18/122
     #88+t1*(t27+t37+t34+t119+t120+t42-t61-t121+2*t122)*t22*t18/1024+(t2
     #7+t37+2*t129+t34+t119+t120+t42-t61-t121)*t1*t7*t18/1024-t1*t149*t2
     #/(-1+l8)*t6*t10*t12*t15*t16/2048-t1*t21*t7*t18/512
      t1 = -1+z
      t2 = -1+l8
      t3 = 1/t2
      t6 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t7 = 1/t6
      t8 = 1/l3*t7
      t10 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t11 = 1/t10
      t12 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t13 = 1/t12
      t15 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t16 = 1/t15
      t17 = f1(1.D0,z,1.D0)
      t19 = t11*t13*t16*t17
      t22 = 1/l2
      t23 = t1*t22
      t29 = dlog(1-l4)
      t30 = dlog(l4)
      t31 = dlog(l1)
      t32 = dlog(2.D0)
      t33 = 2*t32
      t35 = dlog(-t1)
      t36 = 4*t35
      t37 = f3(1.D0,z)
      t39 = dlog(1-t37)
      t40 = f2(1.D0,z)
      t41 = dlog(t40)
      t42 = dlog(z)
      t43 = 2*t42
      t44 = dlog(l2)
      t52 = dlog(l3)
      t59 = 0.3141592653589793D1**2
      t62 = t32**2
      t64 = t42**2
      t94 = -t59-6*t29*t41+12*t62+12*t64-6*t30*t41-12*t41*t32-6*t31*t41+
     #12*t41*t42+6*t31*t39+6*t29*t39+6*t30*t39-48*t35*t42+48*t35*t32+24*
     #t30*t35-12*t30*t42-12*t31*t42+24*t31*t35+6*t29*t30
      t97 = t41**2
      t99 = t39**2
      t117 = t35**2
      t121 = t29**2
      t123 = t31**2
      t127 = t30**2
      t133 = -12*t29*t42+3*t97+3*t99+6*t30*t31+12*t30*t32+12*t29*t32+24*
     #t29*t35+6*t29*t31+12*t31*t32-24*t42*t32-12*t39*t42+48*t117-6*t39*t
     #41+3*t121+3*t123+12*t39*t32+3*t127+24*t35*t39-24*t35*t41
      t150 = dlog(-t2)
      tncoff2 = t1*t3*t8*t19/1024+t23*t3*t7*t19/1024+t1*(t29+t30+t31+t33
     #+t36+t39-t41-t43+2*t44)*t22*t7*t19/1024+(t29+t30+2*t52+t31+t33+t36
     #+t39-t41-t43)*t1*t8*t19/1024-t1*(t94+t133)*t7*t19/12288-t1*(t39*l8
     #-t41*l8+2*t32*l8+t29*l8-2*t42*l8+t31*l8+4*t35*l8+t150*l8+t30*l8)/l
     #8*t3*t7*t11*t13*t16*t17/2048-t23*t8*t19/512
      t1 = -1+z
      t2 = -1+l8
      t3 = 1/t2
      t5 = 1/l8
      t9 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t11 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t14 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t16 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t19 = f1(1.D0,z,1.D0)
      t21 = 1/t9/t11/t14/t16*t19
      t30 = dlog(l1)
      t32 = dlog(-t1)
      t34 = f3(1.D0,z)
      t36 = dlog(1-t34)
      t38 = dlog(1-l4)
      t40 = dlog(-t2)
      t41 = dlog(l4)
      t42 = dlog(l8)
      t43 = f2(1.D0,z)
      t44 = dlog(t43)
      t45 = dlog(2.D0)
      t47 = dlog(z)
      tncoff3 = -t1*t3*t5/l3*t21/1024-t1*t5/l2*t3*t21/1024-t1*(-t30-4*t3
     #2-t36-t38-t40-t41-t42+t44-2*t45+2*t47)*t5*t3*t21/2048
      t1 = -1+z
      t4 = -2+l3
      t5 = 1/t4
      t6 = 1/l3
      t7 = t5*t6
      t10 = 1-l3/2
      t11 = Dsquar311(z,l1,0.D0,t10,l4)
      t13 = s2vr2(z,1.D0,1.D0,t10,l4)
      t16 = s23r2(z,1.D0,1.D0,t10,l4)
      t18 = s1vr2(z,1.D0,1.D0,t10,l4)
      t21 = f1(t10,z,1.D0)
      t23 = 1/t11/t13/t16/t18*t21
      t33 = dlog(1-l4)
      t34 = dlog(l4)
      t35 = dlog(l3)
      t38 = dlog(-t4)
      t41 = dlog(-t1)
      t43 = dlog(l1)
      t44 = f2(t10,z)
      t45 = dlog(t44)
      t46 = dlog(z)
      t48 = f3(t10,z)
      t50 = dlog(1-t48)
      tncoff4 = -t1/l8*t7*t23/512-t1*t5/l2*t6*t23/256+t1*(t33+t34+2*t35+
     #2*t38+4*t41+t43-t45-2*t46+t50)*t7*t23/512
      t1 = -1+z
      t2 = -2+l3
      t3 = 1/t2
      t4 = t1*t3
      t5 = 1/l3
      t11 = 1-l3/2
      t12 = Dsquar311(z,l1,0.D0,t11,l4)
      t14 = s2vr2(z,1.D0,1.D0,t11,l4)
      t17 = s23r2(z,1.D0,1.D0,t11,l4)
      t19 = s1vr2(z,1.D0,1.D0,t11,l4)
      t22 = f1(t11,z,1.D0)
      t24 = 1/t12/t14/t17/t19*t22
      t33 = dlog(1-l4)
      t34 = dlog(l4)
      t35 = dlog(l3)
      t38 = dlog(-t2)
      t41 = dlog(-t1)
      t43 = dlog(l1)
      t44 = f2(t11,z)
      t45 = dlog(t44)
      t46 = dlog(z)
      t48 = f3(t11,z)
      t50 = dlog(1-t48)
      tncoff5 = t4*t5/(-1+l8)*t24/512-t4/l2*t5*t24/256+t1*(t33+t34+2*t35
     #+2*t38+4*t41+t43-t45-2*t46+t50)*t3*t5*t24/512
      t13 = 1-l3/2
      t14 = Dsquar311(z,l1,0.D0,t13,l4)
      t16 = s2vr2(z,1.D0,1.D0,t13,l4)
      t19 = s23r2(z,1.D0,1.D0,t13,l4)
      t21 = s1vr2(z,1.D0,1.D0,t13,l4)
      t24 = f1(t13,z,1.D0)
      tncoff6 = -(-1+z)/l8/(-2+l3)/l3/(-1+l8)/t14/t16/t19/t21*t24/512
      t1 = -1+z
      t4 = 1/l2
      t5 = Dsquar311(z,l1,l2,1.D0,l4)
      t6 = 1/t5
      t7 = t4*t6
      t10 = 1-l1*l2
      t11 = 1-l2
      t12 = s2vr2(z,t10,t11,1.D0,l4)
      t14 = s23r2(z,t10,t11,1.D0,l4)
      t17 = s1vr2(z,t10,t11,1.D0,l4)
      t19 = f1(1.D0,z,t10)
      t21 = 1/t12/t14/t17*t19
      t25 = dlog(1-l4)
      t26 = dlog(l4)
      t27 = dlog(l1)
      t28 = dlog(t11)
      t30 = dlog(-t1)
      t32 = dlog(t10)
      t33 = dlog(2.D0)
      t35 = f2(1.D0,z)
      t36 = dlog(t35)
      t37 = dlog(z)
      t39 = f3(1.D0,z)
      t43 = dlog(1-t39+t39*l1*l2)
      t44 = dlog(l2)
      t46 = t25+t26+t27+t28+4*t30+t32+2*t33-t36-2*t37+t43+2*t44
      tncoff7 = t1/l8*t7*t21/1024-t1*t46*t7*t21/1024+t1*t4/l3*t6*t21/512
      t1 = -1+z
      t2 = 1/l2
      t3 = t1*t2
      t6 = Dsquar311(z,l1,l2,1.D0,l4)
      t7 = 1/t6
      t11 = 1-l1*l2
      t12 = 1-l2
      t13 = s2vr2(z,t11,t12,1.D0,l4)
      t15 = s23r2(z,t11,t12,1.D0,l4)
      t18 = s1vr2(z,t11,t12,1.D0,l4)
      t20 = f1(1.D0,z,t11)
      t22 = 1/t13/t15/t18*t20
      t26 = dlog(1-l4)
      t27 = dlog(l4)
      t28 = dlog(l1)
      t29 = dlog(t12)
      t31 = dlog(-t1)
      t33 = dlog(t11)
      t34 = dlog(2.D0)
      t36 = f2(1.D0,z)
      t37 = dlog(t36)
      t38 = dlog(z)
      t40 = f3(1.D0,z)
      t44 = dlog(1-t40+t40*l1*l2)
      t45 = dlog(l2)
      t47 = t26+t27+t28+t29+4*t31+t33+2*t34-t37-2*t38+t44+2*t45
      tncoff8 = -t3/(-1+l8)*t7*t22/1024-t1*t47*t2*t7*t22/1024+t3/l3*t7*t
     #22/512
      t9 = Dsquar311(z,l1,l2,1.D0,l4)
      t12 = 1-l1*l2
      t13 = 1-l2
      t14 = s2vr2(z,t12,t13,1.D0,l4)
      t17 = s23r2(z,t12,t13,1.D0,l4)
      t19 = s1vr2(z,t12,t13,1.D0,l4)
      t22 = f1(1.D0,z,t12)
      tncoff9 = (-1+z)/l8/l2/(-1+l8)/t9/t14/t17/t19*t22/1024
      t10 = 1-l3/2
      t11 = Dsquar311(z,l1,l2,t10,l4)
      t14 = 1-l1*l2
      t15 = 1-l2
      t16 = s2vr2(z,t14,t15,t10,l4)
      t19 = s23r2(z,t14,t15,t10,l4)
      t21 = s1vr2(z,t14,t15,t10,l4)
      t24 = f1(t10,z,t14)
      tncoff10 = (-1+z)/(-2+l3)/l2/l3/t11/t16/t19/t21*t24/256
      t10 = 1-l3/2
      t11 = Dsquar311(z,l1,l2,t10,l4)
      t14 = 1-l1*l2
      t15 = 1-l2
      t16 = s2vr2(z,t14,t15,t10,l4)
      t19 = s23r2(z,t14,t15,t10,l4)
      t21 = s1vr2(z,t14,t15,t10,l4)
      t24 = f1(t10,z,t14)
      tncoff11 = (-1+z)/(-2+l3)/l2/l3/t11/t16/t19/t21*t24/256
      t1 = -1+z
      t2 = f3(1.D0,z)
      t4 = dlog(1-t2)
      t6 = dlog(l4)
      t8 = dlog(l1)
      t11 = dlog(1-l4)
      t14 = dlog(-t1)
      t16 = f2(1.D0,z)
      t17 = dlog(t16)
      t19 = dlog(z)
      t21 = dlog(2.D0)
      t25 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t26 = 1/t25
      t28 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t30 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t33 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t35 = f1(1.D0,z,1.D0)
      t37 = 1/t28/t30/t33*t35
      tncoff12 = -t1*(-6*t4-6*t6-6*t8-6*t11-24*t14+6*t17+12*t19-12*t21)*
     #t26*t37/12288-t1/l2*t26*t37/1024-t1/l3*t26*t37/1024-t1/l8*t26*t37/
     #2048
      t1 = -1+z
      t4 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t5 = 1/t4
      t7 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t9 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t12 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t14 = f1(1.D0,z,1.D0)
      t16 = 1/t7/t9/t12*t14
      t24 = f3(1.D0,z)
      t26 = dlog(1-t24)
      t28 = dlog(l4)
      t30 = dlog(l1)
      t33 = dlog(1-l4)
      t36 = dlog(-t1)
      t38 = f2(1.D0,z)
      t39 = dlog(t38)
      t41 = dlog(z)
      t43 = dlog(2.D0)
      tncoff13 = -t1/l2*t5*t16/1024-t1/l3*t5*t16/1024-t1*(-6*t26-6*t28-6
     #*t30-6*t33-24*t36+6*t39+12*t41-12*t43)*t5*t16/12288+t1/(-1+l8)*t5*
     #t16/2048
      t6 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t10 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t12 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t15 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t17 = f1(1.D0,z,1.D0)
      tncoff14 = -(-1+z)/l8/(-1+l8)/t6/t10/t12/t15*t17/2048
      t7 = 1-l3/2
      t8 = Dsquar311(z,l1,0.D0,t7,l4)
      t12 = s2vr2(z,1.D0,1.D0,t7,l4)
      t14 = s23r2(z,1.D0,1.D0,t7,l4)
      t17 = s1vr2(z,1.D0,1.D0,t7,l4)
      t19 = f1(t7,z,1.D0)
      tncoff15 = -(-1+z)/(-2+l3)/l3/t8/t12/t14/t17*t19/512
      t7 = 1-l3/2
      t8 = Dsquar311(z,l1,0.D0,t7,l4)
      t12 = s2vr2(z,1.D0,1.D0,t7,l4)
      t14 = s23r2(z,1.D0,1.D0,t7,l4)
      t17 = s1vr2(z,1.D0,1.D0,t7,l4)
      t19 = f1(t7,z,1.D0)
      tncoff16 = -(-1+z)/(-2+l3)/l3/t8/t12/t14/t17*t19/512
      t4 = Dsquar311(z,l1,l2,1.D0,l4)
      t8 = 1-l1*l2
      t9 = 1-l2
      t10 = s2vr2(z,t8,t9,1.D0,l4)
      t12 = s23r2(z,t8,t9,1.D0,l4)
      t15 = s1vr2(z,t8,t9,1.D0,l4)
      t17 = f1(1.D0,z,t8)
      tncoff17 = (-1+z)/l2/t4/t10/t12/t15*t17/1024
      t4 = Dsquar311(z,l1,l2,1.D0,l4)
      t8 = 1-l1*l2
      t9 = 1-l2
      t10 = s2vr2(z,t8,t9,1.D0,l4)
      t12 = s23r2(z,t8,t9,1.D0,l4)
      t15 = s1vr2(z,t8,t9,1.D0,l4)
      t17 = f1(1.D0,z,t8)
      tncoff18 = (-1+z)/l2/t4/t10/t12/t15*t17/1024
      t2 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t10 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t13 = f1(1.D0,z,1.D0)
      tncoff19 = -(-1+z)/t2/t5/t8/t10*t13/2048
      t2 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t10 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t13 = f1(1.D0,z,1.D0)
      tncoff20 = -(-1+z)/t2/t5/t8/t10*t13/2048

CCC Numerator functions

      tcnum1 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,l8,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = CDnum1e0(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = CDnum1e0(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = CDnum1e0(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,l8,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff9)
      tcnum10 = CDnum1e0(var,ndim,z,ux,1-l1*l2,1-l2,1-l3/2,l4,l5,l6,l7,0
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff10)
      tcnum11 = CDnum1e0(var,ndim,z,ux,1-l1*l2,1-l2,1-l3/2,l4,l5,l6,l7,1
     #.D0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff11)
      tcnum12 = CDnum1e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff12)
      tcnum13 = CDnum1e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff13)
      tcnum14 = CDnum1e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-99
     #9.D0,-999.D0,-999.D0,Tcoeff*tncoff14)
      tcnum15 = CDnum1e1(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff15)
      tcnum16 = CDnum1e1(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff16)
      tcnum17 = CDnum1e1(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,0.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff17)
      tcnum18 = CDnum1e1(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,1.D
     #0,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff18)
      tcnum19 = CDnum1e2(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,0.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff19)
      tcnum20 = CDnum1e2(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,1.D0,-
     #999.D0,-999.D0,-999.D0,Tcoeff*tncoff20)

CCC Sum

      CD2e0int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9+tcnum10*tncoff10+tcnum11*tncoff11+tcnum12*tncoff12+t
     #cnum13*tncoff13+tcnum14*tncoff14+tcnum15*tncoff15+tcnum16*tncoff16
     #+tcnum17*tncoff17+tcnum18*tncoff18+tcnum19*tncoff19+tcnum20*tncoff
     #20

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function CD2e1int(var,ndim,Tcoeff)

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
      double precision CDnum1e0
      external CDnum1e0
      double precision CDnum1e1
      external CDnum1e1
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l5 = var(7)
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = var(9)
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t1 = -1+z
      t2 = dlog(l1)
      t3 = dlog(2.D0)
      t6 = dlog(-t1)
      t8 = f3(1.D0,z)
      t10 = dlog(1-t8)
      t11 = f2(1.D0,z)
      t12 = dlog(t11)
      t13 = dlog(z)
      t16 = dlog(1-l4)
      t17 = dlog(l4)
      t20 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t21 = 1/t20
      t23 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t25 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t28 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t30 = f1(1.D0,z,1.D0)
      t32 = 1/t23/t25/t28*t30
      tncoff1 = t1*(t2+2*t3+4*t6+t10-t12-2*t13+t16+t17)*t21*t32/2048-t1/
     #l8*t21*t32/2048-t1/l3*t21*t32/1024-t1/l2*t21*t32/1024
      t1 = -1+z
      t2 = dlog(l1)
      t3 = dlog(2.D0)
      t6 = dlog(-t1)
      t8 = f3(1.D0,z)
      t10 = dlog(1-t8)
      t11 = f2(1.D0,z)
      t12 = dlog(t11)
      t13 = dlog(z)
      t16 = dlog(1-l4)
      t17 = dlog(l4)
      t20 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t21 = 1/t20
      t23 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t25 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t28 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t30 = f1(1.D0,z,1.D0)
      t32 = 1/t23/t25/t28*t30
      tncoff2 = t1*(t2+2*t3+4*t6+t10-t12-2*t13+t16+t17)*t21*t32/2048+t1/
     #(-1+l8)*t21*t32/2048-t1/l3*t21*t32/1024-t1/l2*t21*t32/1024
      t6 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t10 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t12 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t15 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t17 = f1(1.D0,z,1.D0)
      tncoff3 = -(-1+z)/l8/(-1+l8)/t6/t10/t12/t15*t17/2048
      t7 = 1-l3/2
      t8 = Dsquar311(z,l1,0.D0,t7,l4)
      t12 = s2vr2(z,1.D0,1.D0,t7,l4)
      t14 = s23r2(z,1.D0,1.D0,t7,l4)
      t17 = s1vr2(z,1.D0,1.D0,t7,l4)
      t19 = f1(t7,z,1.D0)
      tncoff4 = -(-1+z)/(-2+l3)/l3/t8/t12/t14/t17*t19/512
      t7 = 1-l3/2
      t8 = Dsquar311(z,l1,0.D0,t7,l4)
      t12 = s2vr2(z,1.D0,1.D0,t7,l4)
      t14 = s23r2(z,1.D0,1.D0,t7,l4)
      t17 = s1vr2(z,1.D0,1.D0,t7,l4)
      t19 = f1(t7,z,1.D0)
      tncoff5 = -(-1+z)/(-2+l3)/l3/t8/t12/t14/t17*t19/512
      t4 = Dsquar311(z,l1,l2,1.D0,l4)
      t8 = 1-l1*l2
      t9 = 1-l2
      t10 = s2vr2(z,t8,t9,1.D0,l4)
      t12 = s23r2(z,t8,t9,1.D0,l4)
      t15 = s1vr2(z,t8,t9,1.D0,l4)
      t17 = f1(1.D0,z,t8)
      tncoff6 = (-1+z)/l2/t4/t10/t12/t15*t17/1024
      t4 = Dsquar311(z,l1,l2,1.D0,l4)
      t8 = 1-l1*l2
      t9 = 1-l2
      t10 = s2vr2(z,t8,t9,1.D0,l4)
      t12 = s23r2(z,t8,t9,1.D0,l4)
      t15 = s1vr2(z,t8,t9,1.D0,l4)
      t17 = f1(1.D0,z,t8)
      tncoff7 = (-1+z)/l2/t4/t10/t12/t15*t17/1024
      t2 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t10 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t13 = f1(1.D0,z,1.D0)
      tncoff8 = -(-1+z)/t2/t5/t8/t10*t13/2048
      t2 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t10 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t13 = f1(1.D0,z,1.D0)
      tncoff9 = -(-1+z)/t2/t5/t8/t10*t13/2048

CCC Numerator functions

      tcnum1 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)
      tcnum3 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,l8,-999
     #.D0,-999.D0,-999.D0,Tcoeff*tncoff3)
      tcnum4 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,0.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff4)
      tcnum5 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,1-l3/2,l4,l5,l6,l7,1.D0,
     #-999.D0,-999.D0,-999.D0,Tcoeff*tncoff5)
      tcnum6 = CDnum1e0(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,0.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff6)
      tcnum7 = CDnum1e0(var,ndim,z,ux,1-l1*l2,1-l2,1.D0,l4,l5,l6,l7,1.D0
     #,-999.D0,-999.D0,-999.D0,Tcoeff*tncoff7)
      tcnum8 = CDnum1e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff8)
      tcnum9 = CDnum1e1(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff9)

CCC Sum

      CD2e1int = tcnum1*tncoff1+tcnum2*tncoff2+tcnum3*tncoff3+tcnum4*tnc
     #off4+tcnum5*tncoff5+tcnum6*tncoff6+tcnum7*tncoff7+tcnum8*tncoff8+t
     #cnum9*tncoff9

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      double precision function CD2e2int(var,ndim,Tcoeff)

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
      double precision CDnum1e0
      external CDnum1e0
      NF = 5.0d0
      Pi = 3.141592654d0
      z3 = 1.202056903d0
      buf = 0.00000000001d0
      jac = 1.0d0-2.0d0*buf
      l1 = jac*var(1)+buf
      l2 = jac*var(2)+buf
      l3 = jac*var(3)+buf
      l4 = (dsin(Pi*var(4)/2.0d0))**2
      l5 = var(7)
      l6 = (dsin(Pi*var(8)/2.0d0))**2
      l7 = var(9)
      l8 = jac*var(10)+buf
      z = jac*var(5)+buf
      ux = jac*var(6)+buf
      u = (l3+z-l3*z)/(1-l3+l3*z)
      r = dsqrt(u)
      t = dsqrt(z)
      secsum = 0.0d0


CCC Coefficients

      t2 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t10 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t13 = f1(1.D0,z,1.D0)
      tncoff1 = -(-1+z)/t2/t5/t8/t10*t13/2048
      t2 = Dsquar311(z,l1,0.D0,1.D0,l4)
      t5 = s2vr2(z,1.D0,1.D0,1.D0,l4)
      t8 = s23r2(z,1.D0,1.D0,1.D0,l4)
      t10 = s1vr2(z,1.D0,1.D0,1.D0,l4)
      t13 = f1(1.D0,z,1.D0)
      tncoff2 = -(-1+z)/t2/t5/t8/t10*t13/2048

CCC Numerator functions

      tcnum1 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,0.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff1)
      tcnum2 = CDnum1e0(var,ndim,z,ux,1.D0,1.D0,1.D0,l4,l5,l6,l7,1.D0,-9
     #99.D0,-999.D0,-999.D0,Tcoeff*tncoff2)

CCC Sum

      CD2e2int = tcnum1*tncoff1+tcnum2*tncoff2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

